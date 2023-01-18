struct LocalOp2D{Axis,T}
    data::LocalOperator{T}
end

function capture_operator_or_constant(x)
    if @capture(x, X_[y_Symbol])
        σ = paulidict[X]
        return :(LocalOp2D{$y}($σ, 0:0))
    elseif @capture(x, X_[y_ + i_] | X_[y_ i_])
        σ = paulidict[X]
        return :(LocalOp2D{$y}($σ, ($i):($i)))
    elseif @capture(x, X_[y_ - i_] | X_[y_ -i_])
        σ = paulidict[X]
        return :(LocalOp2D{$y}($σ, (-$i):(-$i)))
    elseif @capture(x, *(A_Symbol, B__))
        if A === :im
            return x
        else
            return :(*($(esc(A)), $(B...)))
        end
    elseif @capture(x, -X_)
        return :(-1 * $X)
    else
        return x
    end
end

function parse_term(ex)
    return MacroTools.postwalk(x -> capture_operator_or_constant(x), ex)
end

function parse_expression(ex)
    terms = []
    if @capture(ex, +(exprs__)) # First, check for a summation of terms
        # Deal with subtractions within summations
        for exp in exprs
            if @capture(exp, X_ - Y_)
                # Convert to a sum
                push!(terms, X, :(-1 * $Y))
            else
                push!(terms, exp)
            end
        end
    elseif @capture(ex, X_ - Y_) # Check for a top level subtraction
        push!(terms, X, :(-1 * $Y))
    else # Otherwise assume a single term
        push!(terms, ex)
    end

    terms = parse_term.(terms) # Parse each term 

    # Initialise vectors to hold terms on each axis
    x_terms = []
    y_terms = []

    # For each term, get it's axis, replace with generic LocalOperator and add to the
    # appropriate vector
    for term in terms
        term_axis = Symbol()
        count = 1
        cleaned_term = MacroTools.postwalk(term) do x
            # Any instance of LocalOp2D{axis} needs to be replaced.
            if @capture(x, LocalOp2D{axis_})
                # Check if all axis labels are the same within a term, but ignore first axis
                # label.
                if !(count == 1) && !(axis === term_axis)
                    # If not, throw an error.
                    throw(
                        ArgumentError(
                            "Differing axis labels $term_axis and $axis used in same term",
                        ),
                    )
                end

                # Store axis information for later use.
                term_axis = axis

                count = count + 1

                # Remove any information about the axis so we can use LocalOperator
                # operations
                return :(LocalOperator)
            else
                # Do nothing on everything else
                return x
            end
        end
        if term_axis === :x
            push!(x_terms, cleaned_term)
        elseif term_axis === :y
            push!(y_terms, cleaned_term)
        else
            push!(x_terms, cleaned_term)
            push!(y_terms, cleaned_term)
        end
    end
    # Add a zero matrix to everything incase of empty vectors.
    # null = LocalOperator(zeros(ComplexF64, 2, 2), 0)
    # return :((+($(x_terms...)))), :((+($(y_terms...))))
    return x_terms, y_terms
end

function hamiltonian(ex)
    terms_x, terms_y = parse_expression(ex)

    # If there exists no y terms then treat x as a normal index
    if length(terms_x) == 0
        terms_x = terms_y
    elseif length(terms_y) == 0
        terms_y = terms_x
    end
    # println(terms_y)
    return quote begin
            rv_x = +($(terms_x...))
            rv_y = +($(terms_y...))
            tuple(rv_x, rv_y)
        end
    end
end

macro hamiltonian(ex)
    return hamiltonian(ex)
end

macro liouvillian(ex_H, ex_Ls...)
    return liouvillian(ex_H, ex_Ls...)
end

function vec2sum(terms)
    return :(+($(terms...)))
end

function liouvillian(ex_H, ex_Ls...)
    terms_H_x, terms_H_y = parse_expression(ex_H)
    ex_H_x = vec2sum(terms_H_x)
    ex_H_y = vec2sum(terms_H_y)

    # tuple of vectors of terms
    ex_Ls_xy = parse_expression.(ex_Ls)

    ex_Ls_x = :([$(vec2sum.(getindex.(ex_Ls_xy, 1))...)])
    ex_Ls_y = :([$(vec2sum.(getindex.(ex_Ls_xy, 2))...)])

    rv = quote
        begin
            r_x = LocalOperators.promote_support($ex_H_x, $(ex_Ls_x)...)
            r_y = LocalOperators.promote_support($ex_H_y, $(ex_Ls_y)...)

            v_x = vectorise(
                LocalOperators.fit($ex_H_x, r_x),
                LocalOperators.fit.($ex_Ls_x, Ref(r_x))...,
            )
            v_y = vectorise(
                LocalOperators.fit($ex_H_y, r_y),
                LocalOperators.fit.($ex_Ls_y, Ref(r_y))...,
            )
            (v_x, v_y)
        end
    end
    return rv
end

macro infliouvillian(ex_H, ex_Ls...)
    H_x, H_y = parse_expression(ex_H)
    ex_H_x = :([$(H_x...)])
    ex_H_y = :([$(H_y...)])
    Ls_xy = parse_expression.(ex_Ls)

    Ls_x = :([$(vec2sum.(getindex.(Ls_xy, 1))...)])
    Ls_y = :([$(vec2sum.(getindex.(Ls_xy, 2))...)])
    rv = quote
        begin
            v_x = infvectorise($ex_H_x, $Ls_x)
            v_y = infvectorise($ex_H_y, $Ls_y)
            (v_x,v_y)
        end
    end
    return rv
end
