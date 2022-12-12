rm_h(ex) = Val(ex.head), ex.args
rm_a1(ex) = Val(ex.args[1]), ex.args[2:end]

@testset "Index parsing" begin
    @test parse_index_pm(Val(:+), :y, 1) == (:y, 1)
    @test parse_index_pm(Val(:-), :i, 3) == (:i, -3)

    @test parse_index(:x) == (:x, 0)
    @test parse_index(:(j + 1)) == (:j, 1)
    @test_throws AssertionError parse_index(:(X[x]))
end

@testset "Single term parsing" begin
    @testset "Positive" begin
        ex0 = :(X[j - 4])
        @test parse_term_single(rm_h(ex0)...) == (σ1, :j, -4)
        ex1 = :(X[i] - X[i])
        @test parse_term_single(rm_h(ex1)...) == (zero(σ0), :i, 0)
        ex2 = :(X[i + 1] - X[i + 1])
        @test parse_term_single(rm_h(ex2)...) == (zero(σ0), :i, 1)
        ex3 = :(X[x - 1] - X[x - 1])
        @test parse_term_single(rm_h(ex3)...) == (zero(σ0), :x, -1)
        ex4 = :(X[i] + Y[i])
        @test parse_term_single(rm_h(ex4)...) == ([0.0 1.0-im; 1.0+im 0.0], :i, 0)
        ex5 = :(Z[i] + Z[i] + Z[i])
        @test parse_term_single(rm_h(ex5)...) == (3 * σ3, :i, 0)
        ex6 = :(Z[i] - Z[i] + Z[i])
        @test parse_term_single(rm_h(ex6)...) == (σ3, :i, 0)
    end
    @testset "Negative" begin
        ex_es = (
            :(X[x] - Y[y]),
            :(X[i] - X[i + 1]),
            :(Z[j] + Y[j] + X[i]),
            :(Z[j] + Y[j] + Y[j - 1]),
        )
        for ex in ex_es
            @test_throws ArgumentError parse_term_single(rm_h(ex)...)
        end
    end
end

@testset "Changing of matrix size" begin
    @test equalise_size([], 0:2) == zeros(ComplexF64, 8, 8)
end
