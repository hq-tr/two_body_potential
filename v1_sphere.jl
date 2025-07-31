module PseudoPotential

using SparseArrays
using Combinatorics
using LinearAlgebra
using Arpack
using WignerSymbols

include("/home/trung/_qhe-julia/Misc.jl")
using .MiscRoutine

function getbasis(filewf::String, N_o::Int64, N_e::Int64)

	ind = 0:N_o-1
	basis_ind = []
	basis = BitVector[]
	dim = 0
	s = (N_o - 1)/2
	m_list = collect(-s:1:s)

	file = open(filewf, "w")
	write(file, "dim" * '\n')
	### Compute all possible basis ###
	for j in combinations(ind, N_e)
		basis_ind = collect(j)
		if sum(m_list[i+1] for i in j) == 0
			basis_dex = dex2bin(basis_ind, N_o)
    		push!(basis, basis_dex)
    		write(file, join(basis_ind) * '\n')
    		write(file, "0.0" * '\n')
    		dim += 1
    	end
    end
    close(file)
#    display(basis)
	return basis, dim
end

function pp_matrix(s::Float64, m::Int)
	dim = Int(2s+1)
	mat = zeros(dim,dim)
	for i in 1:dim
		for j in i:dim
			#println((i,j))
			if abs(i+j-2-2s) < (2s-m)
				mat[i,j] = clebschgordan(Float64,s,i-1-s,s,j-1-s,2s-m,i+j-2-2s)
			end
		end
	end
	return mat
end

function pp_matrix(s::Float64, m::Int32)
	dim = Int(2s+1)
	mat = zeros(dim,dim)
	for i in 1:dim
		for j in i:dim
			#println((i,j))
			if abs(i+j-2-2s) <= (2s-m)
				mat[i,j] = clebschgordan(Float64,s,i-1-s,s,j-1-s,2s-m,i+j-2-2s)
			end
		end
	end
	return mat
end

# This is the new version, calculate the CG coefficients in advance and store in v_mat
function update_element!(H_matrix::SparseMatrixCSC{Float64}, N_o::Int64, i::Int64, j::Int64, basis1::BitVector, basis2::BitVector, v_mat::Matrix{Float64})#; quiet=true)
	if i == j
		basis = findall(basis1)
		H_matrix[i,j] = sum(map(k->abs2(v_mat[k[1],k[2]]), combinations(basis,2)))
	else
		b = basis1 .⊻ basis2
		
		if sum(b) == 4
			m1m2 = findall(basis1 .& b)
			m3m4 = findall(basis2 .& b)
			if sum(m1m2) == sum(m3m4)
				c = count(basis1[m1m2[1]:m1m2[2]])
				d = count(basis2[m3m4[1]:m3m4[2]])
				H_matrix[i,j] = (-1)^(c+d) * v_mat[m1m2[1], m1m2[2]] * v_mat[m3m4[1],m3m4[2]]
				H_matrix[j,i] = H_matrix[i,j]
			end
		end
	end
	return
end


# This is the old version: calculate the CG coefficients on the fly.
function update_element!(H_matrix::SparseMatrixCSC{Float64}, N_o::Int64, i::Int64, j::Int64, basis1::BitVector, basis2::BitVector)#; quiet=true)
	s = (N_o - 1)/2
	b = basis1 .⊻ basis2
	
	if sum(b) == 4
		m1m2 = bin2dex(basis1 .& b)
		m3m4 = bin2dex(basis2 .& b)
		if sum(m1m2) == sum(m3m4)
			#basis1_dex = bin2dex(basis1)
			#basis2_dex = bin2dex(basis2)
			m1 = m1m2[1] - s
			m2 = m1m2[2] - s
			m3 = m3m4[1] - s
			m4 = m3m4[2] - s
			c = count(basis1[m1m2[1]+1:m1m2[2]+1])
			d = count(basis2[m3m4[1]+1:m3m4[2]+1])
			#c = findfirst(x -> x == m2+s, basis1_dex) - findfirst(x -> x == m1+s, basis1_dex) - 1
			#d = findfirst(x -> x == m4+s, basis2_dex) - findfirst(x -> x == m3+s, basis2_dex) - 1
#			println("$s $(m1) $s $(m2) $(2*s-1) $(m1+m2) $s $(m3) $s $(m4) $(2*s-1) $(m3+m4)")
			H_matrix[i,j] = (-1)^(c+d)*clebschgordan(Float64, s, m1, s, m2, (2*s-1), m1+m2) * clebschgordan(Float64, s, m3, s, m4, (2*s-1), m3+m4)
			#H_matrix[i,j] = (-1)^(c+d) * v_mat[m1m2[1], m1m2[2]] * v_mat[m3m4[1],m3m4[2]]
			H_matrix[j,i] = H_matrix[i,j]
		end
	elseif sum(b) == 0
		basis = bin2dex(basis1) .- s
		H_matrix[i,j] = sum(map(k->abs2(clebschgordan(Float64, s, k[1], s, k[2], (2*s-1), k[1]+k[2])), combinations(basis,2)))
	end
	return
end

function v1(N_o::Int64, basis::Vector{BitVector})#; quiet=false)
	dim = length(basis)
	println("The dimension is $(dim)")
	s = (N_o-1)/2
	println("s = $s")
	vmat = pp_matrix(s,1)
	H_matrix = spzeros(dim, dim)
	for i in 1:dim
		#print("\rRow $(i+1)\t\t")
		for j in i:dim
			update_element!(H_matrix, N_o, i, j, basis[i], basis[j],vmat)#; quiet=quiet)
			#print("\r$i\t$j\t\t")
		end
	end
#	display(H_matrix)
	return H_matrix
end

# Below is the full two-body PP code ============================================================

# This is the new version, calculate the CG coefficients in advance and store in v_mat
# Here vmat is a vector of vector corresponds to V_1, V_3, etc...
function update_element!(H_matrix::SparseMatrixCSC{Float64}, N_o::Int64, 
			i::Int64, j::Int64, basis1::BitVector, basis2::BitVector, 
			v_mat::Vector{Matrix{Float64}})
	if i == j
		basis = findall(basis1)
		for v in v_mat
			H_matrix[i,j] += sum(k->abs2(v[k[1],k[2]]), combinations(basis,2))
		end
	else
		b = basis1 .⊻ basis2
		
		if sum(b) == 4
			m1m2 = findall(basis1 .& b)
			m3m4 = findall(basis2 .& b)
			if sum(m1m2) == sum(m3m4)
				c = count(basis1[m1m2[1]:m1m2[2]])
				d = count(basis2[m3m4[1]:m3m4[2]])
				for v in v_mat
					H_matrix[i,j] += (-1)^(c+d) * v[m1m2[1], m1m2[2]] * v[m3m4[1],m3m4[2]]
				end
				H_matrix[j,i] = H_matrix[i,j]
			end
		end
	end
	return
end

# This is the matrix for two-body interaction given a basis
function two_body(N_o::Int64, basis::Vector{BitVector},
				v_list::Vector{Int32}, c_list::Vector{Float64})
	dim = length(basis)
	println("The dimension is $(dim)")
	s = (N_o-1)/2
	println("s = $s")
	vmat = [pp_matrix(s,v_list[i]) * sqrt(c_list[i]) for i in 1:length(v_list)]
	H_matrix = spzeros(dim, dim)
	for i in 1:dim
		#print("\rRow $(i+1)\t\t")
		for j in i:dim
			update_element!(H_matrix, N_o, i, j, basis[i], basis[j],vmat)#; quiet=quiet)
			#print("\r$i\t$j\t\t")
		end
	end
#	display(H_matrix)
	return H_matrix
end

c₊(s::Number,m::Number) = √(s*(s+1) - m*(m+1))
c₋(s::Number,m::Number) = √(s*(s+1) - m*(m-1))
# This is the matrix for L^+ L^- given a basis
function L⁺L⁻(basis::Vector{BitVector})
	dim = length(basis)
	N_o = length(basis[1])
	println("The dimension is $(dim)")
	s = (N_o-1)/2
	println("s = $s")
	H_matrix = spzeros(dim, dim)
	for i in 1:dim
		#print("\rRow $(i+1)\t\t")
		# Diagonal term]
		print("\r$i\t$i\t\t")
		possibles = findall(basis[i][2:end] .& .!basis[i][1:end-1])
		H_matrix[i,i] = sum(m->c₋(s,m-s)c₊(s,m-1-s), possibles)
		for j in (i+1):dim
			print("\r$i\t$j\t\t")
			b = basis[i] .⊻ basis[j]
			
			if sum(b) == 4
				ms = findall(b) .- 1
				if abs(ms[2]-ms[1])==1
					H_matrix[i,j] = c₊(s,ms[1]-s)*c₋(s,ms[end]-s)
					H_matrix[j,i] = H_matrix[i,j]
				end
			end
		end
	end
#	display(H_matrix)
	return H_matrix
end

LplusLminus = L⁺L⁻
LL = L⁺L⁻


# The following two functions do similar things to the two above.
# Except instead of being registered as matrix entries, the terms are added to a sum
# This sum is the variational energy of a given state.
function update_energy!(ϵ::Vector{Float64}, N_o::Int64, 
			i::Int64, j::Int64, basis1::BitVector, basis2::BitVector, 
			v_mat::Vector{Matrix{Float64}}, coefs::Vector{T} where T <: Number)
	if i == j
		basis = findall(basis1)
		for v in v_mat
			ϵ .+= conj(coefs[i]) * coefs[j] * sum(map(k->abs2(v[k[1],k[2]]), combinations(basis,2)))
		end
	else
		b = basis1 .⊻ basis2
		
		if sum(b) == 4
			m1m2 = findall(basis1 .& b)
			m3m4 = findall(basis2 .& b)
			if sum(m1m2) == sum(m3m4)
				c = count(basis1[m1m2[1]:m1m2[2]])
				d = count(basis2[m3m4[1]:m3m4[2]])
				for v in v_mat
					ϵ .+= 2*real(conj(coefs[i]) * coefs[j] * (-1)^(c+d) * v[m1m2[1], m1m2[2]] * v[m3m4[1],m3m4[2]])
				end
			end
		end
	end
	#println(ϵ)
	return
end

# There are two functions for evaluating the variational energy
# They implement slightly different methods
function eval_energy(N_o::Int64, 
			i::Int64, j::Int64, basis1::BitVector, basis2::BitVector, 
			v_mat::Vector{Matrix{Float64}}, coef::Number)
	ϵ = 0
	if i == j
		basis = findall(basis1)
		for v in v_mat
			ϵ += coef * sum(k->abs2(v[k[1],k[2]]), combinations(basis,2))
		end
	else
		b = basis1 .⊻ basis2
		
		if sum(b) == 4
			m1m2 = findall(basis1 .& b)
			m3m4 = findall(basis2 .& b)
			if sum(m1m2) == sum(m3m4)
				c = count(basis1[m1m2[1]:m1m2[2]])
				d = count(basis2[m3m4[1]:m3m4[2]])
				for v in v_mat
					ϵ += 2*real(coef * (-1)^(c+d) * v[m1m2[1], m1m2[2]] * v[m3m4[1],m3m4[2]])
				end
			end
		end
	end
	return ϵ
end

function two_body_energy(N_o::Int64, basis::Vector{BitVector}, 
				coefs::Vector{T} where T <: Number,
				v_list::Vector{Int32}, c_list::Vector{Float64})
	dim = length(basis)
	println("The dimension is $(dim)")
	s = (N_o-1)/2
	println("s = $s")
	println("Calculating the CG coefficients in advance")
	energy = [0.0]
	# The single-particle pseudopotential matrix is modified by the square root of the coefficient.
	# This is because the energy contribution is always the product of two entries.
	@time vmat = [pp_matrix(s,v_list[i]).*sqrt(c_list[i]) for i in 1:length(v_list)]
	#energy = sum(k->eval_energy(N_o, k[1],k[2], basis[k[1]], basis[k[2]], vmat, conj(coefs[k[1]])*coefs[k[2]]), with_replacement_combinations(1:dim,2))
	for k in with_replacement_combinations(1:dim,2)
		update_energy!(energy, N_o, k[1],k[2], basis[k[1]], basis[k[2]], vmat, coefs)
	end
	#display(H_matrix)
	return energy[1]
end


export v1, two_body, two_body_energy, L⁺L⁻, LplusLminus, LL
end
