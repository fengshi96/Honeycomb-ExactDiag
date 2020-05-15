using LinearAlgebra
include("Lattice.jl")

function Hamiltonian()
	Kzz = 1; Kyy = 1; Kxx = 1 ##########################change later
	nsite, mesh_, nn_, indx_, indy_= Honeycomb(2,2,"PBC")
	KxxConn_ = fill!(Matrix{Float64}(undef, nsite, nsite),0)
	KyyConn_ = fill!(Matrix{Float64}(undef, nsite, nsite),0)
	KzzConn_ = fill!(Matrix{Float64}(undef, nsite, nsite),0)
	
	for i::Int8 in 1:nsite
		# Kxx_Conn
		j::Int8 = nn_[i,1]
		if i < j && j>=1
			KxxConn_[i,j] = Kxx
			KxxConn_[j,i] = Kxx
		end
	
		# Kyy_Conn
		j = nn_[i,2]
		if i < j && j>=1
			KyyConn_[i,j] = Kyy
			KyyConn_[j,i] = Kyy
		end	

		# Kzz_Conn
		j = nn_[i,3]
		if i < j && j>=1
			KzzConn_[i,j] = Kzz
			KzzConn_[j,i] = Kzz
		end
	
	end
	
	




























end

Hamiltonian()
#show(stdout, "text/plain", nn_); println()

