
macro liouvillian(ex_H, ex_Ls...)
    H_x, Ls_x, H_y, Ls_y = parse_lindblad(ex_H, ex_Ls...)

    L_x = vectorise(H_x, Ls_x...)
    L_y = vectorise(H_y, Ls_y...)

    return L_x::Matrix{ComplexF64}, L_y::Matrix{ComplexF64}
end

parse_pauli(::Val{:X}) = pauli[1]
parse_pauli(::Val{:Y}) = pauli[2]
parse_pauli(::Val{:Z}) = pauli[3]
parse_pauli(::Val{:I}) = pauli[4]
parse_pauli(::Val{:P}) = σp
parse_pauli(::Val{:M}) = σm

function parse_lindblad(ex_H, ex_Ls...)
    dict_H = parse_lin(ex_H)
    dict_Ls = parse_lin.(ex_Ls)

    terms_H_x, terms_Ls_x = equalise_sizes(dict_H[:x], getindex.(dict_Ls, :x)...)
    terms_H_y, terms_Ls_y = equalise_sizes(dict_H[:y], getindex.(dict_Ls, :y)...)
    return sum(terms_H_x), sum.(terms_Ls_x), sum(terms_H_y), sum.(terms_Ls_y)
end

function get_indexrange(vecs)
    mins = get_rangemin.(vecs)
    maxs = get_rangemax.(vecs)

    min_i = min(mins...)
    max_i = max(maxs...)

    return Int(min_i):Int(max_i)
end

function get_rangemin(vec)
    if length(vec) == 0
        return Inf
    end
    min_i = min((@. first(Tuple(keys(vec))))...)
    return min_i
end
function get_rangemax(vec)
    if length(vec) == 0
        return -Inf
    end
    max_i = max((@. last(Tuple(keys(vec))))...)
    return max_i
end

function equalise_sizes(vec_H, vec_Ls...)
    args = (vec_H, vec_Ls...)
    r = get_indexrange(args)
    terms_H = equalise_size(vec_H, r)
    terms_Ls = equalise_size.(vec_Ls, Ref(r))
    return terms_H, terms_Ls
end

function equalise_size(vec, r)
    if length(vec) == 0
        return [zeros(ComplexF64, (2^length(r), 2^length(r)))]
    end

    terms = []

    for dict in vec
        term_i = Vector{Matrix{ComplexF64}}([])
        for i in r
            push!(term_i, dict[i])
        end
        push!(terms, kron(1, term_i...))
    end
    return terms
end

function parse_lin(ex)
    dict = Dict(:x => [], :y => [])
    return parse_expr!(dict, ex, 1)
end

function parse_expr!(dict, ex, pm::Int)
    # Check if expression is :call or :ref 
    if ex.head === :call
        # If :call, then get the operator
        op = ex.args[1]
        if op === :+ # If +, then parse the terms in the sum
            for exp in ex.args[2:end]
                parse_expr!(dict, exp, pm)
            end
        elseif op === :- # If -, then parse the two terms, negating the second
            exp_p = ex.args[2]
            exp_m = ex.args[3]
            parse_expr!(dict, exp_p, pm)
            parse_expr!(dict, exp_m, pm * -1)
        elseif op === :* # If *, then we have reached a term and can append it to the dict
            term = parse_term(Val(:*), ex.args[2:end])
            mult_term!(term, pm)
            push_term!(dict, term)
        end
    elseif ex.head === :ref
        term = parse_term(Val(:ref), ex.args)
        mult_term!(term, pm)
        push_term!(dict, term)
    end
    return dict
end

function push_term!(dict, term)
    if !all(y -> y[2] === term[1][2], term)
        throw(ArgumentError("Different indices used in term"))
    else
        axis = term[1][2]
    end

    term_kv = [sym[3] => sym[1] for sym in term]

    term_dict = DefaultOrderedDict(pauli[4], term_kv...)

    if axis === :x
        push!(dict[:x], term_dict)
    elseif axis === :y
        push!(dict[:y], term_dict)
    else
        push!(dict[:x], term_dict)
        push!(dict[:y], term_dict)
    end

    return dict
end

function mult_term!(term, pm)
    term[1] = (term[1][1] * pm, term[1][2], term[1][3])
    return term
end

function parse_term(::Val{:call}, all_args)
    return parse_term(Val(all_args[1]), all_args[2:end])
end
function parse_term(::Val{:ref}, args)
    return [parse_term_single(Val(:ref), args)]
end
# Returns the expression as a vector of terms
function parse_term(::Union{Val{:+},Val{:-}}, args)
    terms = parse_term.(Val.(getfield.(args, :head)), getfield.(args, :args))
    return terms
end

# Returns a single term as a vector of pauli's
function parse_term(::Val{:*}, args)
    mult = 1
    if typeof(args[1]) <: Number
        mult = args[1]
        start = 2
    elseif args[1] == :im
        mult = im
        start = 2
    else
        start = 1
    end
    terms =
        parse_term_single.(
            Val.(getfield.(args[start:end], :head)), getfield.(args[start:end], :args)
        )
    return mult_term!(terms, mult)
end

# Expression of the form X[x ± i] 
function parse_term_single(::Val{:ref}, args)
    σ = parse_pauli(Val(args[1]))
    xy, i = parse_index(args[2])
    return σ, xy, i
end
# Expression of the form X[x + i] + X[x + i]
function parse_term_single(::Val{:call}, args)
    return parse_term_single(Val(args[1]), args[2:end])
end

function parse_term_single(::Val{:*}, args)
    mult = 1
    if typeof(args[1]) <: Number 
        mult = args[1]
        start = 2
    elseif args[1] == :im
        mult = im
        start = 2
    else
        start = 1
    end
    mults =
        parse_term_single.(
            Val.(getfield.(args[start:end], :head)), getfield.(args[start:end], :args)
        )
    σs = getindex.(mults, 1)
    xys = getindex.(mults, 2)
    pms = getindex.(mults, 3)
    if !(all(y -> y == xys[1], xys))
        throw(ArgumentError("Cannot multiply over different axis within a term"))
    end
    if !(all(y -> y == pms[1], pms))
        throw(ArgumentError("Cannot mutliply over different indices within a term"))
    end
    return *(mult, σs...), xys[1], pms[1]
end

function parse_term_single(::Val{:+}, args)
    sums = parse_term_single.(Val.(getfield.(args, :head)), getfield.(args, :args))
    σs = getindex.(sums, 1)
    xys = getindex.(sums, 2)
    pms = getindex.(sums, 3)
    if !(all(y -> y == xys[1], xys))
        throw(ArgumentError("Cannot sum over different axis within a term"))
    end
    if !(all(y -> y == pms[1], pms))
        throw(ArgumentError("Cannot sum over different indices within a term"))
    end
    return sum(σs), xys[1], pms[1]
end
function parse_term_single(::Val{:-}, args)
    sums = parse_term_single.(Val.(getfield.(args, :head)), getfield.(args, :args))
    σs = getindex.(sums, 1)
    xys = getindex.(sums, 2)
    pms = getindex.(sums, 3)
    if !(all(y -> y == xys[1], xys))
        throw(ArgumentError("Cannot subtract over different axis within a term"))
    end
    if !(all(y -> y == pms[1], pms))
        throw(ArgumentError("Cannot subtract over different indices within a term"))
    end
    return σs[1] - σs[2], xys[1], pms[1]
end

### Index parsing
function parse_index(ex::Expr)
    @assert ex.head === :call
    args = ex.args
    return parse_index_pm(Val(args[1]), args[2], args[3])
end

parse_index(symb::Symbol) = symb, 0
parse_index_pm(::Val{:+}, xy, pm) = xy, pm
parse_index_pm(::Val{:-}, xy, pm) = xy, -pm
