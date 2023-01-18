# TrotterMacros.jl

Simple macros for parsing infinite, local Hamiltonians and Lindblad equations expressed as sums of Pauli matrices into a (super-)operator matrix representation of the associated trotter gate. 

## Overview

This package exports two convenience macros: `@hamiltonian` and `@liouvillian`. The former constructs a $l$-local trotter gate corresponding to the Hamiltonian supplied, where $l$ is the maximum locality of the terms in the sum. For example, the $2$-local trotter gate corresponding to the transverse-field Ising model can be constructed like
```julia
h_x, _ @hamiltonian Z[i]*Z[i+1] + X[i]
```
The `Z` corresponds to the pauli-$z$ matrix, with the other pauli matrices denoted similarly. The indices are purely relative within each term, that is `Z[i]*Z[i+1]*Z[i+2]` is parsed as $\sigma^z \otimes \sigma^z \otimes \sigma^z$, `X[i-1]*X[i] + Y[i]*Y[i+1]` is parsed as $\sigma^x \otimes \sigma^x  + \sigma^y \otimes \sigma^y$. As such, the overall locality of the resulting trotter gate is dictated by the maximum and the minimum integer appearing across the terms. Whether one starts indexing at $0$, $-1$ or $34$ has no effect on the output to the output. 

Each macro is by default two-dimensional. That is, one can use the special index labels `x` and `y` to specify neighbour interactions along only the respective axis. For example,
```julia
h_x, h_y = @hamiltonian Z[x-1]*X[x]*Z[x+1] + X[y-1]*Z[y]*X[y+1]
```
results in two different trotter gates, one for the $x$-axis, and one for the $y$-axis. Any index other than `x` and `y` will apply to both axes. If one only cares about a single dimension, then the second returned value can be discarded by writing:
```
h, = @hamiltonian ... #or
h, _ = @hamiltonian ...
```
The `@liouvillian` macro works much the same. Any number of Lindblad jump operators can be included as so:
```julia
l_x, l_y = @liouvillian Z[i]*Z[i+1] P[i] M[i] # or
l_x, l_y = @liouvillian(Z[i]*Z[i+1], P[i], M[i])
```
The output is then the *vectorised* form of the local super-operator corresponding to the equationsupplied to the macro. As Julia utilises column-major ordering, we vectorise according to $A x B \rightarrow (B^{T} \otimes A ) \vec{x}$. 
