using LinearAlgebra

function Honeycomb(LLX, LLY, BC::String = "OBC")

	#LLX = 2; LLY = 2
	Nsite = LLX * LLY * 2
	number1neigh = 3
	println("creating Honeycomb lattice... ")

	scalex::Float64 = 2; scaley::Float64 = 4.0/sqrt(3)
	t1 = [1.0 * scalex; 0]
	t2 = [0.5 * scalex; sqrt(3)/2.0 * scaley]

	indx_ = Vector{Int8}(undef, Nsite)  # x coordinate of sites
	indy_ = Vector{Int8}(undef, Nsite)  # y coordinate of sites
	mesh_ = fill!(Matrix{Int8}(undef, LLX*2+LLY,LLY*2),0)  # lattice on meshgrid
	nn_ = Matrix{Int8}(undef, Nsite, number1neigh)  # nearest neighbors
	
	# ---------------- Construct Lattice mesh ------------------------ 
	# ---------------- Construct Lattice mesh ------------------------ 
	# ---------------- Construct Lattice mesh ------------------------ 
	xv::Int8 = 1  # xv is unit cell length in x (vertival) direction
	counter::Int8 = 1	
	for i in 1: LLX
		if i!=1
			xv += t1[1]
		end
				
		xa::Int8 = xv; xb::Int8 = xv + 1
		ya::Int8 = 1; yb::Int8 = 2
		
		for j in 1:LLY
			indx_[counter] = xa
			indy_[counter] = ya
			mesh_[xa,ya] = counter
			#println("xa=", xa, " ya=",ya, " counter=",counter)
			counter += 1
			
			indx_[counter] = xb
			indy_[counter] = yb
			mesh_[xb,yb] = counter
			#println("xb=", xb, " yb=",yb, " counter=",counter)
			counter += 1
			
			xa = xv + j*t2[1]; xb = (xv + 1) + j*t2[1]
			ya = 1 + j*t2[2]; yb = 2 + j*t2[2]
										
		end		
	end

	# ------------ Construct matrix of 1st n.n.s --------------------- 
	# ------------ Construct matrix of 1st n.n.s --------------------- 
	# ------------ Construct matrix of 1st n.n.s --------------------- 
	println("Looking for nearest neighbors... ")
	xmax = findmax(indx_)[1]   # findmax returns both max value and its index
	ymax = findmax(indy_)[1]
	
	for i in 1:Nsite
		ix::Int8 = indx_[i]  # coordinate of n-th site in matrix
		iy::Int8 = indy_[i]
				
		#----------------------------OBC-----------------------------------
		#  n.n in x-bond
		jx::Int8 = ix + 1; jy::Int8 = iy + 1  # move 1 step forward in x = (1,1) direction
		if jx <= xmax && jy <= ymax && mesh_[jx,jy] != 0
			j = mesh_[jx,jy]  # site index of n.n. in x direction
			nn_[i,1] = j
			nn_[j,1] = i 	
		end
	
		#  n.n in y-bond
		jx = ix + 1; jy = iy - 1  # move 1 step in y = (1,-1) direction		
		if jx <= xmax && jy <= ymax && jy >= 1 && mesh_[jx,jy] != 0
			j = mesh_[jx,jy]  # site index of n.n. in x direction
			nn_[i,2] = j
			nn_[j,2] = i 	
		end	
				
		#  n.n in z-bond
		jx = ix; jy = iy + 1  # move 1 step in z = (0,1) direction		
		if jx <= xmax && jy <= ymax && mesh_[jx,jy] != 0
			j = mesh_[jx,jy]  # site index of n.n. in x direction
			nn_[i,3] = j
			nn_[j,3] = i 	
		end
		#----------------------------OBC-----------------------------------
		
		
		#--------------------------Apply PBC-------------------------------	
		if BC == "PBC"
			# z-bond
			jx = ix - LLY 
			jy = 1
			if jx >= 1 && iy == ymax && mesh_[jx,jy]!=0
				j = mesh_[jx,jy]
				nn_[i,3] = j
				nn_[j,3] = i
			end
			
			# y-bond
			jx = ix + 2 * LLX - 1 
			jy = iy + 1
			if jx <= xmax && iy <= ymax && iy%2 !=0 && mesh_[jx,jy]!=0
				j = mesh_[jx,jy]
				nn_[i,2] = j
				nn_[j,2] = i
			end
		end			
	
	end
	
	return mesh_, nn_

end

mesh_, nn_ = Honeycomb(2,2,"PBC")	
#show(stdout, "text/plain", nn_); println()
#show(stdout, "text/plain", mesh_); println()

