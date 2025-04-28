include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("v1_sphere.jl")
using .FQH_states
using .PseudoPotential
using LinearAlgebra
using Arpack
using SparseArrays
using ArgMacros

using BenchmarkTools

using Random

function main()
    @inlinearguments begin
        @argumentrequired String fname "-f" "--file-name-root"
        @argumentrequired Int dim "-d" "--dimension"
        @argumentrequired String intname "-i" "--interaction-file"
        @argumentflag noeigenstate "--no-eigenstate"
    end

    println("============")
    if dim < 1
        println("The dimension must be positive. Terminating.")
        return false
    else
        println("IMPORTANT: all the input states must contain the same monomial basis listed in the same order.")
        println("The input files must be named in the following format: <name root><index>")
        println("where <index> is an integer from 0 to dim-1.")

        if isfile(fname * "0")
            state = readwf(fname * "0")
            basis = state.basis
            h_dim = length(basis)
            coefs = zeros(Float64,h_dim,dim)
        else
            println("File name '$(fname)0' not found. Terminating.")
            return false
        end

        for i in 1:dim
            if isfile(fname * string(i-1))
                state = readwf(fname * string(i-1))
                coefs[:,i] = identity.(state.coef)
            else
                print("File name '$(fname)$(i-1)' not found. Terminating.")
                return false
            end
        end

        N_o = length(basis[1])
        println("$(N_o) orbitals.")
    end


    #@time basis, dim = getbasis(filewf, N_o, N_e)

    v_list = Int32[]
    c_list = Float64[]

    if isfile(intname)
        open(intname) do f
            for line in map(s->split(s),readlines(f))
                append!(v_list,parse(Int32,line[1]))
                append!(c_list,parse(Float64,line[2]))
            end
        end
    else
        print("Interaction file '$(intname)' not found. Terminating.")
        return false
    end

    println("============")
    println("Constructing the Hamiltonian")

    #@time ϵ = two_body_energy(N_o, basis, coefs, v_list, c_list)

    @time H_matrix = two_body(N_o, basis, v_list, c_list)
    display(H_matrix)

    println("--------")

    if dim == 1
        @time ϵ = coefs' * H_matrix * coefs
        println("The energy is $ϵ")
        println("--------")
    else
        @time mat = coefs' * H_matrix * coefs
        ϵ = eigvals(mat)
        println("The eigenvalues are")
        display(ϵ)
    end

    # Save the eigenstates
    if !noeigenstate
        for i in 1:dim
            state_coef = coefs * vecs[:,i]
            state = FQH_state(basis, state_coef)
            printwf(state;"g_$(i-1)")
            i+=1
        end
    end

    return true
#    println(transpose(ϕ))

end

@time main()