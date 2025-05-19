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
        @argumentrequired String fname "-f" "--file-name"
        @argumentrequired String intname "-i" "--interaction-file"
    end
    println("Calculate two-body pseudopotential variational energy of a given state.")

    state = readwf(fname)

    basis = state.basis
    coefs = state.coef

    N_o = length(basis[1])
    println("$(N_o) orbitals.")

    #@time basis, dim = getbasis(filewf, N_o, N_e)


    v_list = Int32[]
    c_list = Float64[]

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

    println("--------")
    println("Constructing the Hamiltonian")

    @time ϵ = two_body_energy(N_o, basis, coefs, v_list, c_list)

    #@time H_matrix = two_body(N_o, basis, v_list, c_list)

    #println("--------")
    #println(size(coefs))
    #println(size(H_matrix))

    #@time ϵ = coefs' * H_matrix * coefs



    println("The energy is $ϵ")
    println("--------")


#    println(transpose(ϕ))

end

@time main()