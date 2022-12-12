# TrotterMacros.jl

Simple macros for parsing Hamiltonians and Lindblad equations expressed as sums of Pauli matrices into a (super-)operator matrix representation of the associated trotter gate. 

## Assumptions

The macros can parse sums and subtractions of terms of form:
$
   \lambda_i \Sigma_{i \pm l_1} \Sigma_{i \pm l_2} \cdot \Sigma_{i \pm l_k} 
$
where
$
    \Sigma  = 
$
That is, for each term within the summation, factor out all constants and factor out any summations *within* this term that contain sums over different indices. 
