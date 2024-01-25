include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("v1_sphere.jl")
using .FQH_states
using .PseudoPotential
using LinearAlgebra
using Arpack
using SparseArrays

using BenchmarkTools

using Random

function main()
    println("Calculate two-body pseudopotential variational energy of a given state.")


    println("Input state file name: ")

    state = readwf(readline())

    basis = state.basis
    coefs = state.coef

    N_o = length(basis[1])
    println("$(N_o) orbitals.")

    #@time basis, dim = getbasis(filewf, N_o, N_e)

    println("Input m for Vₘ and the corresponding coefficient. ")
    println("Each pp term takes one line, with two numbers separated by a space.")
    println("Put a 0 to end")


    v_list = Int32[]
    c_list = Float64[]

    reading = true
    while reading
        data = readline()
        if data == "0"
            reading = false
        else
            try
                pp = split(data)
                push!(v_list,parse(Int32, pp[1]))
                push!(c_list,parse(Float64,pp[2]))
            catch
                println("Invalid input. Try again or input 0 to end.")
            end
        end
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