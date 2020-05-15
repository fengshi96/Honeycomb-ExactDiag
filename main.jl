using SparseArrays
using LinearAlgebra
using Arpack
include("src/Hamiltonian.jl")

LLX = 2; LLY = 2
@time Ham = Hamiltonian(LLX,LLY,"PBC")
@time eigvals, eigvecs = eigs(Ham, nev = 2, which=:SR) # SR = Smallest Real

println("\nEigen values:"); show(stdout, "text/plain", real.(eigvals)); println()





