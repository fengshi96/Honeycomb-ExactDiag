"""
A LatticeBond is a struct which represents
a single bond in a geometrical lattice or
else on interaction graph defining a physical
model such as a quantum Hamiltonian.
LatticeBond has the following data fields:
* s1::Int -- number of site 1
* s2::Int -- number of site 2
* x1::Float64 -- x coordinate of site 1
* y1::Float64 -- y coordinate of site 1
* x2::Float64 -- x coordinate of site 2
* y2::Float64 -- y coordinate of site 2
* type::String -- optional description of bond type
"""
struct LatticeBond
  s1::Int
  s2::Int
  x1::Float64
  y1::Float64
  x2::Float64
  y2::Float64
  type::String
end

"""
    LatticeBond(s1::Int,s2::Int)
    LatticeBond(s1::Int,s2::Int,
                x1::Real,y1::Real,
                x2::Real,y2::Real,
                type::String="")
Construct a LatticeBond struct by
specifying just the numbers of sites
1 and 2, or additional details including
the (x,y) coordinates of the two sites and
an optional type string.
"""
function LatticeBond(s1::Int, s2::Int)
  return LatticeBond(s1, s2, 0.0, 0.0, 0.0, 0.0, "")
end

function LatticeBond(s1::Int, s2::Int, bondtype::String="")
  return LatticeBond(s1, s2, 0.0, 0.0, 0.0, 0.0, bondtype)
end

function LatticeBond(
  s1::Int, s2::Int, x1::Real, y1::Real, x2::Real, y2::Real, bondtype::String=""
)
  cf(x) = convert(Float64, x)
  return LatticeBond(s1, s2, cf(x1), cf(y1), cf(x2), cf(y2), bondtype)
end

"""
Lattice is an alias for Vector{LatticeBond}
"""
const Lattice = Vector{LatticeBond}

"""
    triangular_lattice(Nx::Int,
                       Ny::Int;
                       kwargs...)::Lattice
Return a Lattice (array of LatticeBond
objects) corresponding to the two-dimensional
triangular lattice of dimensions (Nx,Ny).
By default the lattice has open boundaries,
but can be made periodic in the y direction
by specifying the keyword argument
`yperiodic=true`. The lattice is enumerated 
in y-major 
"""
function triangular_lattice(Nx::Int, Ny::Int; kwargs...)::Tuple{Lattice, Matrix{Int64}}
  yperiodic = get(kwargs, :yperiodic, false)
  yperiodic = yperiodic && (Ny > 2)
  N = Nx * Ny
  Nbond = 3N - 2Ny + (yperiodic ? 0 : -2Nx + 1)
  latt = Lattice(undef, Nbond)
  geometry = Matrix{Int64}(undef, (N, 3))  # columbs are: index, x, y
  b = 0
  for n in 1:N
    x = div(n - 1, Ny) + 1  # x coordinate in square graph
    y = mod(n - 1, Ny) + 1  # y coordinate

    geometry[n, 1] = n
    geometry[n, 2] = x
    geometry[n, 3] = y

    # x-direction bonds
    if x < Nx
      latt[b += 1] = LatticeBond(n, n + Ny, "a")
    end

    # 2d bonds
    if Ny > 1
      # vertical / y-periodic diagonal bond
      if (n + 1 <= N) && ((y < Ny) || yperiodic)
        latt[b += 1] = LatticeBond(n, n + 1, "b")
      end
      # periodic vertical bond
      if yperiodic && y == 1
        latt[b += 1] = LatticeBond(n, n + Ny - 1, "b")
      end
      # diagonal bonds
      if x < Nx && y < Ny
        latt[b += 1] = LatticeBond(n, n + Ny + 1, "c")
      end
    end
  end
  return latt, geometry
end
