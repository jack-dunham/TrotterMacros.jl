coherent(H::AbstractMatrix) = -im * (kron(one(H), H) - kron(transpose(H), one(H)))

function dissipator(L::AbstractMatrix)

    rv =
        1 / 2 * (
            2 * conj(L) ⊗ L - one(L) ⊗ (L' * L) -
            (transpose(L) * conj(L)) ⊗ one(L)
        )

    return rv
end

function vectorise(H, Ls...)
    D = sum(dissipator.(Ls))
    return coherent(H) + D
end
