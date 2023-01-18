using TrotterMacros
import TrotterMacros: hamiltonian, liouvillian, pauli, coherent, dissipator, generate_terms
using Test
using TestSetExtensions

@testset ExtendedTestSet "All tests" begin
    σx, σy, σz, σ0 = pauli
    @testset "Hamiltonians (expression parser)" begin
        @testset "Single term" begin
            @testset "Single site index" begin
                x, y = @hamiltonian Y[i]
                @test x == y
                @test size(x) == (2, 2)
                @test x == σy

                x, y = @hamiltonian Y[i] * Y[i]
                @test x == y
                @test size(x) == (2, 2)
                @test x == σ0

                x, y = @hamiltonian Y[i + 1] * Y[i + 1]
                @test x == y
                @test size(x) == (2, 2)
                @test x == σ0
            end

            @testset "Index syntax" begin
                x, y = @hamiltonian Z[i] * Z[i + 1]
                @test x == y
                @test size(x) == (4, 4)
                @test x == kron(σz, σz)

                # Hcat style
                x, y = @hamiltonian Z[i] * Z[i +1]
                @test x == y
                @test size(x) == (4, 4)
                @test x == kron(σz, σz)
                x, y = @hamiltonian Z[i] * Z[i -1]
                @test x == y
                @test size(x) == (4, 4)
                @test x == kron(σz, σz)
                # Negative indices
                x, y = @hamiltonian Z[i - 1] * X[i] * Z[i + 1]
                @test x == y
                @test size(x) == (8, 8)
                @test x == kron(σz, σx, σz)
            end
            @testset "Constant multiplication" begin
                # Constant multiplication
                x, y = @hamiltonian 3 * Z[i - 1] * X[i] * Z[i + 1]
                @test x == y
                @test size(x) == (8, 8)
                @test x == 3 * kron(σz, σx, σz)

                # Imaginary multiplication
                x, y = @hamiltonian im * Z[i - 1] * X[i] * Z[i + 1]
                @test x == y
                @test size(x) == (8, 8)
                @test x == im * kron(σz, σx, σz)

                # Interpolation
                J = 3
                x, y = @hamiltonian J * Z[i - 1] * X[i] * Z[i + 1]
                @test x == y
                @test size(x) == (8, 8)
                @test x == J * kron(σz, σx, σz)
            end
            @testset "Sums within terms" begin
                x, y = @hamiltonian Z[i - 1] * (X[i] + Y[i] + Z[i]) * Z[i + 1]
                @test x == y
                @test size(x) == (8, 8)
                @test x == kron(σz, σx + σy + σz, σz)
                x, y = @hamiltonian Z[i - 1] * (X[i] - Y[i]) * Z[i + 1]
                @test x == y
                @test size(x) == (8, 8)
                @test x == kron(σz, σx - σy, σz)
                x, y = @hamiltonian Z[i - 1] * (-X[i] - Y[i]) * Z[i + 1]
                @test x == y
                @test size(x) == (8, 8)
                @test x == kron(σz, -σx - σy, σz)
                x, y = @hamiltonian Z[i - 1] * (X[i] - 2 * Y[i]) * Z[i + 1]
                @test x == y
                @test size(x) == (8, 8)
                @test x == kron(σz, σx - 2 * σy, σz)
                x, y = @hamiltonian Z[i - 1] * (X[i] - Y[i + 2]) * Z[i + 1]
                @test x == y
                @test size(x) == (16, 16)
                @test x == kron(σz, σx, σz, σ0) - kron(σz, σ0, σz, σy)
            end
            @testset "Single axis edge cases" begin
                # Single axis
                x, y = @hamiltonian 3 * Z[x - 1] * X[x] * Z[x + 1]
                @test x == 3 * kron(σz, σx, σz)
                @test y == x
            end
        end
        @testset "Two terms" begin
            # Top level addition
            x, y = @hamiltonian X[i] * X[i + 1] + Y[i] * Y[i + 1]
            @test x == y
            @test size(x) == (4, 4)
            @test x == kron(σx, σx) + kron(σy, σy)

            # Top level subtraction
            x, y = @hamiltonian X[i] * X[i + 1] - Y[i] * Y[i + 1]
            @test x == y
            @test size(x) == (4, 4)
            @test x == kron(σx, σx) - kron(σy, σy)

            # Differing axis
            x, y = @hamiltonian X[x] * X[x + 1] - X[y] * X[y + 1]
            @test x == -y
            @test size(x) == size(y) == (4, 4)
            @test x == kron(σx, σx) == -y
        end
    end
    @testset "Liouvillians (vectorising)" begin
        @testset "Utils" begin
            H = rand(ComplexF64, 4, 4)
            L = rand(ComplexF64, 2, 2)
            rv = generate_terms(TrotterMacros.LocalOperator(H, 1:2), 2)
            @test length(rv) == 1
            @test rv[1] ≈ H

            rv2 = generate_terms(TrotterMacros.LocalOperator(H, 1:2), 4)
            @test length(rv2) == 3
            @test rv2[1] ≈ kron(H, one(H), one(H))
            @test rv2[2] ≈ kron(one(H), H, one(H))
            @test rv2[3] ≈ kron(one(H), one(H), H)
        end
        @testset "Vectorising" begin
            H = rand(ComplexF64, 2, 2)
            L = rand(ComplexF64, 2, 2)
            p = rand(ComplexF64, 2, 2)
            Hp = -im * (H * p - p * H)
            D = L * p * L' - 0.5 * (L' * L * p + p * L' * L)

            Hv = coherent(H)
            pv = reshape(p, 4)

            @test size(Hv) == (4, 4)
            @test reshape(Hv * pv, (2, 2)) ≈ Hp
            @test reshape(dissipator(L) * pv, (2, 2)) ≈ D
            @test coherent(kron(H, one(H)) + kron(one(H), H)) ≈
                sum(coherent.([kron(H, one(H)), kron(one(H), H)]))

            H1 = rand(ComplexF64, 2, 2)
            H2 = rand(ComplexF64, 2, 2)
            p1 = rand(ComplexF64, 2, 2)
            p2 = rand(ComplexF64, 2, 2)

            H = kron(H1, H2)
            p = kron(p1, p2)
            pv = reshape(p, 16)

            Hp = -im * (H * p - p * H)
            Hv = coherent(H)

            @test reshape(Hv * pv, (4, 4)) ≈ Hp
        end
        @testset "Dissipative ising mode" begin
            H = kron(σz, σz)
            L = 1 / 2 * (σx - im * σy) #\sigma-
            L2 = kron(L, one(L))

            Hv = coherent(H)
            Dv = dissipator(L2)

            x, y = @liouvillian Z[i] * Z[i + 1] M[i]
            @test x == Hv + Dv

            x, y = @infliouvillian Z[i] * Z[i + 1] M[i]
            xn =
                coherent(kron(σz, σz)) +
                1 / 2 * dissipator(kron(L, one(L))) +
                1 / 2 * dissipator(kron(one(L), L))

            @test x ≈ xn

            x, y = @liouvillian Z[i] * Z[i + 1]
            @test x == coherent(kron(σz, σz))
        end
    end
end
