using TrotterMacros
import TrotterMacros: pauli, σ0, σ1, σ2, σ3, σp, σm
import TrotterMacros: equalise_size
import TrotterMacros: parse_index_pm, parse_index, parse_term_single, parse_term
using Test
using TestSetExtensions

@testset ExtendedTestSet "All tests" begin
    @includetests ARGS
end

