export LocalOperator
using LinearAlgebra

struct LocalOperator{T<:Number} <: AbstractMatrix{T}
    data::Matrix{T}
    support::UnitRange{Int}
    localdim::Int
    function LocalOperator(data::AbstractMatrix{T}, support::UnitRange, localdim::Int) where {T}
        D = localdim
        n, m = size(data)
        if !(n == m)
            throw(ArgumentError("matrix size $(size(data)) is not square"))
        elseif !(isinteger(log(D, n)))
            throw(ArgumentError("matrix size $(size(data)) is not compatible with local dimension D=$D"))
        else
            return new{T}(data,support, localdim)
        end
    end
end

const LocalOp = LocalOperator

struct LocalDimensionMismatch <: Exception
    a::LocalOp
    b::LocalOp
end

function Base.showerror(io::IO, e::LocalDimensionMismatch)
    print(io, "LocalDimensionMismatch: local dimensions must match: a has local dim $(localdim(e.a)), b has local dim $(localdim(e.b))")
end
# Default has local dimension 2
LocalOperator(a::AbstractMatrix, r::UnitRange) = LocalOp(a,r, 2)

Base.size(A::LocalOp) = size(A.data)

support(A::LocalOp) = A.support
minsupport(A::LocalOp) = support(A)[begin]
maxsupport(A::LocalOp) = support(A)[end]

localdim(A::LocalOp) = A.localdim

function Base.showarg(io::IO, A::LocalOp, toplevel)
    return print(io, "$(localdim(A))-local ",typeof(A), " on sites ", minsupport(A), " to ", maxsupport(A))
end

Base.getindex(A::LocalOp, i::Int) = getindex(A.data, i)
Base.setindex!(A::LocalOp, v, i::Int) = setindex!(A.data, v, i)

Base.IndexStyle(::Type{<:LocalOp}) = IndexStyle(Matrix{ComplexF64})

Base.similar(A::LocalOp) = LocalOp(similar(A.data), support(A), localdim(A))
function Base.similar(A::LocalOp, ::Type{<:Number}, dims::Dims) 
    D = localdim(A)
    size_change = Int(log(D,dims[1]) - log(D,size(A)[1]))
    r = support(A) .+ size_change
    return LocalOp(similar(A.data, ComplexF64, dims), r, D)
end

# Fallback 
Base.similar(A::LocalOp, dims::Tuple{Base.OneTo{Int64}}) = similar(A.data, dims)

function promote_support(a::LocalOp, b::LocalOp)
    localdim(a) == localdim(b) || throw(LocalDimensionMismatch(a,b))
    l = min(minsupport(a), minsupport(b))
    u = max(maxsupport(a), maxsupport(b))
    return l:u
end

function promote_support(a::LocalOp, b::LocalOp, c::LocalOp, args::Vararg{<:LocalOp})
    return promote_support(promote_support(a, b), c, args...)
end

# Pad `a` on the left and right with `kl` and `kr` 2x2 identity matrices respectively.
function pad(a::LocalOp, kl::Int, kr::Int)
    return kron(1, fill(one(a), kl)..., a, fill(one(a), kr)...)
end

function fit(a::LocalOp, r::UnitRange{Int})
    kl = minsupport(a) - r[begin]
    kr = r[end] - maxsupport(a)
    return LocalOp(pad(a, kl, kr), r)
end

# Basic operations

Base.:*(x::Number, a::LocalOp) = mul!(similar(a), x, a)
Base.:*(a::LocalOp, x::Number) = mul!(similar(a), a, x)

function Base.:*(a::LocalOp, b::LocalOp) 
    r = promote_support(a, b)
    aa = fit(a, r)
    bb = fit(b, r)

    cc = similar(aa)

    mul!(cc,aa,bb)

    return cc
end

function Base.:+(a::LocalOp, b::LocalOp)
    r = promote_support(a, b)
    aa = fit(a, r)
    bb = fit(b, r)

    axpy!(1, aa, bb)

    return bb
end

Base.:-(a::LocalOp, b::LocalOp) = +(a, -1 * b)

function LinearAlgebra.mul!(a::LocalOp, b::LocalOp, x::Number) 
    localdim(a) == localdim(b) || throw(LocalDimensionMismatch(a,b))
    mul!(a.data, b.data, x)
    return a
end
function LinearAlgebra.mul!(a::LocalOp, x::Number, b::LocalOp) 
    localdim(a) == localdim(b) || throw(LocalDimensionMismatch(a,b))
    mul!(a.data,x, b.data)
    return a
end
function LinearAlgebra.axpy!(x::Number, a::LocalOp, b::LocalOp) 
    localdim(a) == localdim(b) || throw(LocalDimensionMismatch(a,b))
    axpy!(x, a.data, b.data)
    return b
end
function LinearAlgebra.axpby!(x::Number, a::LocalOp, y::Number, b::LocalOp) 
    localdim(a) == localdim(b) || throw(LocalDimensionMismatch(a,b))
    axpy!(x, a.data, y, b.data)
    return b
end

function LinearAlgebra.mul!(c::LocalOp,a::LocalOp,b::LocalOp,x::Number,y::Number)
    mul!(c.data,a.data,b.data,x,y)
    return c
end
