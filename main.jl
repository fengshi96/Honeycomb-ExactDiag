using SparseArrays
using LinearAlgebra
using Arpack
include("src/Hamiltonian.jl")
include("src/Parameter.jl")


param = GetParameter("input.inp")
Nstates = param.Nstates
@time Ham = Hamiltonian(param)
@time eigvals, eigvecs = eigs(Ham, nev = Nstates, which=:SR, tol = 0.001) # SR = Smallest Real

println("\nEigen values:"); show(stdout, "text/plain", real.(eigvals)); println()





