module TrotterMacros

using DataStructures
using MacroTools
using LocalOperators

export @hamiltonian, @liouvillian, @infliouvillian

include("utils.jl")
include("parse.jl")

end # module TrotterMacros
