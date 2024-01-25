include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("v1_sphere.jl")
using .FQH_states
using .PseudoPotential
using LinearAlgebra
using Arpack
using SparseArrays

using BenchmarkTools

function main(fname="")
    println("----------")
    println("Full ED of two-body pseudopotential from a given basis file.")
    println("----------")

    if length(fname)==0
        println("Input basis file name: ")
        fname = readline()
    end    

    basis = readwf(fname).basis

    N_o = length(basis[1])
    println("$(N_o) orbitals.")

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

    #@time basis, dim = getbasis(filewf, N_o, N_e)

    println("--------")
    println("Constructing the Hamiltonian")

    @time H_matrix = two_body(N_o, basis, v_list, c_list)#; quiet=true)

    display(H_matrix)

    println("--------")

    println("Diagonalizing with ARPACK")

    @time λ, ϕ = eigs(H_matrix, which=:SM)

    #display(ϕ)

    println("Eigenvalues = ")
    display(real.(λ))

    println("--------")
    gs_coef = ϕ[:,1]
    println(length(gs_coef))
    ground_state = FQH_state(basis, gs_coef)
    printwf(ground_state;fname="g_$(fname)_0")

    println("Saved ground state as g_$(fname)_0.")


#    println(transpose(ϕ))

end

@time main()