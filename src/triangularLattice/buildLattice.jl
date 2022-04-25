include("myLattice.jl")

let
    lattice, geometry = triangular_lattice(4,4, yperiodic=true)
    for l in lattice
        println(l)
    end
    println("\ngeometry"); show(stdout, "text/plain", geometry); println()
end
