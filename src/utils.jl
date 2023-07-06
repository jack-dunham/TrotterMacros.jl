const σ0 = Matrix{ComplexF64}([1 0; 0 1])
const σ1 = Matrix{ComplexF64}([0 1; 1 0])
const σ2 = Matrix{ComplexF64}([0 -im; im 0])
const σ3 = Matrix{ComplexF64}([1 0; 0 -1])

const σp = 1 / 2 * (σ1 + im * σ2) # Pauli+ matrix
const σm = 1 / 2 * (σ1 - im * σ2) # Pauli- matrix

const pauli = (σ1, σ2, σ3, σ0)
const paulidict = Dict(
    :I => pauli[4], :X => pauli[1], :Y => pauli[2], :Z => pauli[3], :P => σp, :M => σm
)

const ⊗ = kron

coherent(H::AbstractMatrix) = -im * (one(H) ⊗ H - transpose(H) ⊗ one(H))

function dissipator(L::AbstractMatrix)
    rv = conj(L) ⊗ L - 0.5 * (one(L) ⊗ (L' * L) + (transpose(L) * conj(L)) ⊗ one(L))

    return rv
end

function vectorise(H, Ls...)
    D = sum(dissipator.(Ls))
    return coherent(H) + D
end

function generate_terms(op::LocalOperator, k::Int)
    mat = op.data

    l = locality(op)

    num_terms = k - l + 1

    term_vec = fill(one(mat), num_terms)
    term_vec[1] = mat

    rv = [kron(one(eltype(mat)), circshift(term_vec, i - 1)...) for i in 1:num_terms]

    return rv
end
function infvectorise(Hs::Vector, Ls::Vector{<:LocalOperator}; is2d=false)
    k = max(locality.(Hs)..., locality.(Ls)...)
    H_terms = reduce(vcat, scalefunc.(coherent, generate_terms.(Hs, k), k; is2d=is2d))
    L_terms = reduce(vcat, scalefunc.(dissipator, generate_terms.(Ls, k), k; is2d=is2d))

    rv = sum(H_terms) + sum(L_terms)

    return rv
end

function scalefunc(f, v::Vector{<:AbstractMatrix}, k; is2d=false)
    num_terms = length(v)
    if !(num_terms == k)
        is2d = false
    end

    coeff = 2^(is2d)

    @debug "" k length(v) coeff
    
    rv = map(x -> (1 / (length(v) * coeff)) * f(x), v)
    return rv
end
