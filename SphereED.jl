module TwoBodyGet
using WignerSymbols
function pp_matrix(s::Float64, m::Int32)
	dim = Int(2s+1)
	mat = zeros(dim,dim)
	for i in 1:dim
		for j in 1:dim
			#println((i,j))
			if abs(i+j-2-2s) <= (2s-m)
				mat[i,j] = clebschgordan(Float64,s,-i+s+1,s,-j+s+1,2s-m,-i-j+2s+2)
			end
		end
	end
    #@show mat
	return mat
end

function get_twobody(n_orb, v_list::Vector{Int32}, c_list::Vector{Float64})
    twoBody = Vector{Tuple{ComplexF64, Tuple{Int64, Int64,Int64,Int64}}}()
    s = (n_orb-1)/2
    vmat = [pp_matrix(s,v_list[i]) * sqrt(c_list[i]) for i in 1:length(v_list)]
    for i = 1:n_orb
        for j = 1:i
            for l = 1:n_orb
                k = i+j-l
                if k>0 && k<n_orb+1
                    h = 0
                    for v in vmat
                        h += conj(v[j,i])*v[l,k]
                    end
                    push!(twoBody,(h,(i,j,l,k)))
                end
            end
        end
    end
    return twoBody
end
export pp_matrix,get_twobody
end

module OneBodyGet
single_particle_state_sphere(θ::Number, φ::Number, S::Number, m::Number) = cos(θ/2)^(S+m) * sin(θ/2)^(S-m) * exp(m*im*φ) / sphere_coef(S,m)

single_particle_state_sphere(θ::Array{T} where T<: Number, φ::Array{T} where T<: Number, S::Number, m::Number) = cos.(θ./2).^(S.+m) .* sin.(θ/2).^(S.-m) .* exp.(m*im.*φ) / sphere_coef(S,m)

# Normalization coefficient on the sphere
sphere_coef(S,m) =   sqfactorial(S-m)/sqfactorial(S+m+1, 2S+1)


sqfactorial(N) = prod(map(sqrt, 1:N)) 
# square root of N!, avoid overflow for large N and more efficient than sqrt(factorial(big(N)))
# (overflow starts at N=21)

sqfactorial(n,N) = prod(map(sqrt, n:N))
function get_oneBody(n_orb::Int64,θ_list::Vector{Float64}, ϕ_list::Vector{Float64}, V_list::Vector{Float64})
    oneBody = Vector{Tuple{ComplexF64, Tuple{Int64, Int64}}}()
    S = (n_orb-1)/2
    for i = 1:length(θ_list)
        θ = θ_list[i]
        ϕ = ϕ_list[i]
        height = V_list[i]
        sfunction = map(m->single_particle_state_sphere(θ,ϕ,S,m), S:-1:-S)
        for j = 1:n_orb
            for k = 1:j
                if k == j
                    push!(oneBody,(conj(sfunction[j])*sfunction[k]*height/2,(j,k)))
                else
                    push!(oneBody,(conj(sfunction[j])*sfunction[k]*height,(j,k)))
                end
            end
        end
    end
    return oneBody
end
using WignerD
function get_oneBody_widebump(n_orb::Int64,θ_list::Vector{Float64}, ϕ_list::Vector{Float64}, pinsize_list::Vector{Int64})
    oneBody = Vector{Tuple{ComplexF64, Tuple{Int64, Int64}}}()
    S = (n_orb-1)/2
    for j = 1:n_orb
        for k = 1:j
            h_jk = 0.0+0.0im
            for i = 1:length(θ_list)
                θ = θ_list[i]
                ϕ = ϕ_list[i]
                W = WignerD.wignerD(S,0,θ,ϕ)
                pinsize = pinsize_list[i]
        
                h_jk += sum([W[Int(S-m0+1),j]*conj(W[Int(S-m0+1),k]) for m0 in S:-1:(S-pinsize+1)])
            end
            if k == j
                push!(oneBody,(h_jk/2,(j,k)))
                #for m0 = S:-1:S-pinsize+1
                #    push!(oneBody,(W[Int(S-m0+1),j]*conj(W[Int(S-m0+1),k])/2,(j,k)))
                #end
            else
                push!(oneBody,(h_jk,(j,k)))
                #for m0 = S:-1:S-pinsize+1
                #    push!(oneBody,(W[Int(S-m0+1),j]*conj(W[Int(S-m0+1),k]),(j,k)))
                #end
            end
        end
    end
    return oneBody
end

function get_oneBody_gaussian(n_orb::Int64,θ_list::Vector{Float64}, ϕ_list::Vector{Float64}, a::Float64)
    gaussian_coefficient(a::Float64,m::Number) = (a^2 / (a^2 + 1))^(m+1) / (2π*a^2)
    gaussian_coefficient_normalized(a::Float64,m::Number) = gaussian_coefficient(a,m) / gaussian_coefficient(a,0) # make sure m=0 coefficient is 1
    oneBody = Vector{Tuple{ComplexF64, Tuple{Int64, Int64}}}()
    S = (n_orb-1)/2
    for j = 1:n_orb
        for k = 1:j
            h_jk = 0.0+0.0im
            for i = 1:length(θ_list)
                θ = θ_list[i]
                ϕ = ϕ_list[i]
                W = WignerD.wignerD(S,0,θ,ϕ)
        
                h_jk += sum([W[Int(S-m0+1),j]*conj(W[Int(S-m0+1),k])*gaussian_coefficient_normalized(a,S-m0) for m0 in S:-1:(-S)])
            end
            if k == j
                push!(oneBody,(h_jk/2,(j,k)))
            else
                push!(oneBody,(h_jk,(j,k)))
            end
        end
    end
    return oneBody
end

export get_oneBody,get_oneBody_widebump,get_oneBody_gaussian
end


module ConstructManybodyMatrix
using SparseArrays
function calcSign(I::Int64, K::Int64)::Int64
    M = I-(I&(2K-1)) #occupied states to the left of (i)
    btwnCnt = count_ones(M) #counts all occupied to the left of (i)
    sign = (-1)^btwnCnt
    return sign
end

#annihilation operator
function C(i::Int64, I::Int64, sign::Int64)::Tuple{Int64, Int64}
    if sign == 0 #if sign is 0, the many body state is anihilated so can just proceed
        return 0, 0
    end
    K = 2^(i-1)  #integer representation of single particle state
    if (I & K) != K # if (i) was already unoccupied, then anihilate many body state
        return 0, 0
    else
        L = I - K  #anihilate (i) in I if (i) is occupied in I
        sign *= calcSign(I, K)
        return L, sign
    end
end

#creation operator
function CDag(i::Int64, I::Int64, sign::Int64)::Tuple{Int64, Int64}
    if sign == 0 #if sign is 0, the many body state is anihilated so can just proceed
        return 0, 0
    end
    K = 2^(i-1)  #integer representation of single particle state
    if I & K == K # if state is already occupied, anihilate many body state
        return 0, 0
    else
        L = I + K #create (i) in I
        sign *= calcSign(I, K)
        return L, sign
    end
end

#two-body part of Hamiltonian
function calcV(twoBody::Vector{Tuple{ComplexF64, Tuple{Int64, Int64,Int64, Int64}}}, States::Vector{Int64})::SparseMatrixCSC{ComplexF64, Int64}
    Rows = [Vector{Int64}() for _ in 1:Threads.nthreads()+1]
    Cols = [Vector{Int64}() for _ in 1:Threads.nthreads()+1]
    Vals = [Vector{ComplexF64}() for _ in 1:Threads.nthreads()+1]
    Threads.@threads for i in eachindex(twoBody)
        (v, (i, j, l, k)) = twoBody[i]
        for IInd in eachindex(States)
            I = States[IInd]
            sign = 1
            J, sign = C(k, I, sign)
            J, sign = C(l, J, sign)
            J, sign = CDag(j, J, sign)
            J, sign = CDag(i, J, sign)
            if sign == 0 #if State is anihilated
                continue
            end
            JInd = searchsortedfirst(States,J)
            if JInd == (length(States)+1) || States[JInd]!=J
                continue
            end
            push!(Rows[Threads.threadid()], JInd)
            push!(Cols[Threads.threadid()], IInd)
            push!(Vals[Threads.threadid()], sign*v) 
        end
    end
    rows=reduce(vcat,Rows)
    cols=reduce(vcat,Cols)
    vals=reduce(vcat,Vals)
    #FREE UP RAM
    Rows=nothing
    Cols=nothing
    Vals=nothing
    GC.gc() #garbage collect
    V = sparse(rows, cols, vals, length(States),length(States))
    rows=nothing
    cols=nothing
    vals=nothing
    GC.gc() 
    return V
end

function calcV!(V::SparseMatrixCSC{ComplexF64, UInt64},twoBody::Vector{Tuple{ComplexF64, Tuple{Int64, Int64,Int64, Int64}}}, States::Vector{Int64},subzone::Int64;num_zones=10)
    Rows = [Vector{UInt64}() for _ in 1:Threads.nthreads()]
    Cols = [Vector{UInt64}() for _ in 1:Threads.nthreads()]
    Vals = [Vector{ComplexF64}() for _ in 1:Threads.nthreads()]
    base_size,remainder = divrem(length(twoBody),num_zones)
    if subzone <= remainder
        begin_index = (base_size+1)*(subzone-1)+1
        end_index = (base_size+1)*subzone
    else
        begin_index = (base_size+1)*remainder + (subzone-remainder-1)*base_size+1
        end_index = (base_size+1)*remainder + (subzone-remainder)*base_size
    end
    Threads.@threads for i in begin_index:end_index
        (v, (i, j, l, k)) = twoBody[i]
        for IInd in eachindex(States)
            I = States[IInd]
            sign = 1
            J, sign = C(k, I, sign)
            J, sign = C(l, J, sign)
            J, sign = CDag(j, J, sign)
            J, sign = CDag(i, J, sign)
            if sign == 0 #if State is anihilated
                continue
            end
            JInd = searchsortedfirst(States,J)
            if JInd == (length(States)+1) || States[JInd]!=J
                continue
            end
            if VERSION >= v"1.12"
                tid = Threads.nthreads() > 1 ? Threads.threadid() - 1 : Threads.threadid()
            else
                tid = Threads.threadid()
            end
            push!(Rows[tid], JInd)
            push!(Cols[tid], IInd)
            push!(Vals[tid], sign*v) 
        end
    end
    rows=reduce(vcat,Rows)
    cols=reduce(vcat,Cols)
    vals=reduce(vcat,Vals)
    #FREE UP RAM
    Rows=nothing
    Cols=nothing
    Vals=nothing
    GC.gc() #garbage collect
    #V += sparse(rows, cols, vals, length(States),length(States))
    for (i,(row,col)) in enumerate(zip(rows,cols))
        V[row,col] += vals[i]
    end
    rows=nothing
    cols=nothing
    vals=nothing
    GC.gc() 
    return
end

    
function calcT(oneBody::Vector{Tuple{ComplexF64, Tuple{Int64, Int64}}}, States::Vector{Int64})::SparseMatrixCSC{ComplexF64, Int64}
    rows = Vector{Int64}()
    cols = Vector{Int64}()
    vals = Vector{ComplexF64}()
    for value in oneBody
        (t, (i, j)) = value
        for IInd in eachindex(States)
            I = States[IInd]
            sign = 1
            J, sign = C(j, I, sign)
            J, sign = CDag(i, J, sign)
            if sign == 0 #if State is anihilated
                continue
            end
            JInd = searchsortedfirst(States,J) #REQUIRES THAT States BE SORTED. see http://www.jlhub.com/julia/manual/en/function/searchsortedfirst
            if JInd == (length(States)+1) || States[JInd]!=J
                continue
            end
            push!(rows, JInd)
            push!(cols, IInd)
            push!(vals, sign*t)
        end
    end
    T = sparse(rows, cols, vals, length(States),length(States))
    return T
end

function calcT!(T::SparseMatrixCSC{ComplexF64,UInt64},oneBody::Vector{Tuple{ComplexF64, Tuple{Int64, Int64}}}, States::Vector{Int64},subzone::Int64,λ::Float64;num_zones=10)
    Rows = [Vector{UInt64}() for _ in 1:Threads.nthreads()]
    Cols = [Vector{UInt64}() for _ in 1:Threads.nthreads()]
    Vals = [Vector{ComplexF64}() for _ in 1:Threads.nthreads()]
    base_size,remainder = divrem(length(oneBody),num_zones)
    if subzone <= remainder
        begin_index = (base_size+1)*(subzone-1)+1
        end_index = (base_size+1)*subzone
    else
        begin_index = (base_size+1)*remainder + (subzone-remainder-1)*base_size+1
        end_index = (base_size+1)*remainder + (subzone-remainder)*base_size
    end 
    Threads.@threads for i in begin_index:end_index
        (t, (i, j)) = oneBody[i]
        if VERSION >= v"1.12"
            tid = Threads.nthreads() > 1 ? Threads.threadid() - 1 : Threads.threadid()
        else
            tid = Threads.threadid()
        end
        for IInd in eachindex(States)
            I = States[IInd]
            sign = 1
            J, sign = C(j, I, sign)
            J, sign = CDag(i, J, sign)
            if sign == 0 #if State is anihilated
                continue
            end
            JInd = searchsortedfirst(States,J) #REQUIRES THAT States BE SORTED. see http://www.jlhub.com/julia/manual/en/function/searchsortedfirst
            if JInd == (length(States)+1) || States[JInd]!=J
                continue
            end
            push!(Rows[tid], JInd)
            push!(Cols[tid], IInd)
            push!(Vals[tid], sign*t) 
        end
    end
    rows=reduce(vcat,Rows)
    cols=reduce(vcat,Cols)
    vals=reduce(vcat,Vals)
    #FREE UP RAM
    Rows=nothing
    Cols=nothing
    Vals=nothing
    GC.gc() #garbage collect
    for (i,(row,col)) in enumerate(zip(rows,cols))
        T[row,col] += λ*vals[i]
        T[col,row] += λ*conj(vals[i])
    end
    rows=nothing
    cols=nothing
    vals=nothing
    GC.gc() 
    return 
end

function calcT!(T::SparseMatrixCSC{ComplexF64,UInt64},oneBody::Vector{Tuple{ComplexF64, Tuple{Int64, Int64}}}, States::Vector{Int64},subzone::Int64,λ::Float64,n_sectors::Integer;num_zones=10)
    Rows = [Vector{UInt64}() for _ in 1:Threads.nthreads()]
    Cols = [Vector{UInt64}() for _ in 1:Threads.nthreads()]
    Vals = [Vector{ComplexF64}() for _ in 1:Threads.nthreads()]
    base_size,remainder = divrem(length(oneBody),num_zones)
    if subzone <= remainder
        begin_index = (base_size+1)*(subzone-1)+1
        end_index = (base_size+1)*subzone
    else
        begin_index = (base_size+1)*remainder + (subzone-remainder-1)*base_size+1
        end_index = (base_size+1)*remainder + (subzone-remainder)*base_size
    end 
    Threads.@threads for i in begin_index:end_index
        (t, (i, j)) = oneBody[i]
        if abs(i-j) % n_sectors != 0
            continue 
        end
        if VERSION >= v"1.12"
            tid = Threads.nthreads() > 1 ? Threads.threadid() - 1 : Threads.threadid()
        else
            tid = Threads.threadid()
        end
        for IInd in eachindex(States)
            I = States[IInd]
            sign = 1
            J, sign = C(j, I, sign)
            J, sign = CDag(i, J, sign)
            if sign == 0 #if State is anihilated
                continue
            end
            JInd = searchsortedfirst(States,J) #REQUIRES THAT States BE SORTED. see http://www.jlhub.com/julia/manual/en/function/searchsortedfirst
            if JInd == (length(States)+1) || States[JInd]!=J
                continue
            end
            push!(Rows[tid], JInd)
            push!(Cols[tid], IInd)
            push!(Vals[tid], sign*t) 
        end
    end
    rows=reduce(vcat,Rows)
    cols=reduce(vcat,Cols)
    vals=reduce(vcat,Vals)
    #FREE UP RAM
    Rows=nothing
    Cols=nothing
    Vals=nothing
    GC.gc() #garbage collect
    for (i,(row,col)) in enumerate(zip(rows,cols))
        T[row,col] += λ*vals[i]
        T[col,row] += λ*conj(vals[i])
    end
    rows=nothing
    cols=nothing
    vals=nothing
    GC.gc() 
    return 
end

export calcV,calcT,CDag,C,calcV!,calcT!
end


