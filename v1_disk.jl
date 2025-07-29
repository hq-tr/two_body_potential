module PseudoPotential

using SparseArrays
using Combinatorics
using LinearAlgebra
using Arpack
using WignerSymbols

include("/home/trung/_qhe-julia/Misc.jl")
using .MiscRoutine

# Below is the full two-body PP code ============================================================

# Calculate the matrix element given m1,m2,m3,m4

# This is the new version, calculate the CG coefficients in advance and store in v_mat
# Here vmat is a vector of vector corresponds to V_1, V_3, etc...
function update_element!(H_matrix::SparseMatrixCSC{Float64}, N_o::Int64, 
			i::Int64, j::Int64, basis1::BitVector, basis2::BitVector, 
			v_list::Vector{Int32},c_list::Vector{Float64})
	if i == j
		basis = findall(basis1)
		for (v,c) in zip(v_list,c_list)
			H_matrix[i,j] += sum(k->h_factor(k[1]-1,k[2]-1,k[1]-1,k[2]-1,v,c), combinations(basis,2))
		end
	else
		b = basis1 .⊻ basis2
		
		if sum(b) == 4
			m1m2 = findall(basis1 .& b)
			m3m4 = findall(basis2 .& b)
			if sum(m1m2) == sum(m3m4)
				p1= count(basis1[m1m2[1]:m1m2[2]])
				p2 = count(basis2[m3m4[1]:m3m4[2]])
				for (v,c) in zip(v_list,c_list)
					H_matrix[i,j] += (-1)^(p1+p2) * h_factor(m1m2[1]-1,m1m2[2]-1,m3m4[1]-1,m3m4[2]-1,v,c)
				end
				H_matrix[j,i] = H_matrix[i,j]
			end
		end
	end
	return
end

function h_factor(m1,m2,m3,m4,v,c)
	factor = (1/√(2))^(m1+m2+m3+m4) * factorial(big(m1+m2-v)) * factorial(v) * sqfactorial(m1)*sqfactorial(m2)*sqfactorial(m3)*sqfactorial(m4)
	sum = 0.
	for i in 0:min(v,m2), j in 0:min(v,m4)
		ii = v - i
		jj = v - j
		if ii ≤ m1 && jj ≤ m3
			sum += (-1)^(i+j) * c * factor / (factorial(i) * factorial(big(m2-i)) * factorial(ii) * factorial(big(m1-ii)) * factorial(j) * factorial(big(m4-j)) * factorial(jj) * factorial(big(m3-jj)))
		end
	end
	return sum
end

# This is the matrix for two-body interaction given a basis
function two_body(N_o::Int64, basis::Vector{BitVector},
				v_list::Vector{Int32}, c_list::Vector{Float64})
	dim = length(basis)
	println("The dimension is $(dim)")

	H_matrix = spzeros(dim, dim)
	for i in 1:dim
		#print("\rRow $(i+1)\t\t")
		for j in i:dim
			update_element!(H_matrix, N_o, i, j, basis[i], basis[j],v_list,c_list)
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
