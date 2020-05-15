using SparseArrays
using LinearAlgebra
using Arpack
include("src/Hamiltonian.jl")

LLX = 3; LLY = 3
Ham = Hamiltonian(LLX,LLY,"PBC")
eigvals, eigvecs = eigs(Ham, nev = 2, which=:SR) # SR = Smallest Real
println(eigvals)





