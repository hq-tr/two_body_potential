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
    println("Full ED of V₁ pseudopotential from a given basis file.")
    println("----------")

    if length(fname)==0
        println("Input basis file name: ")
        fname = readline()
    end    

    basis = readwf(fname).basis

    N_o = length(basis[1])
    println("$(N_o) orbitals.")

    #@time basis, dim = getbasis(filewf, N_o, N_e)

    println("--------")
    println("Constructing the Hamiltonian")

    @time H_matrix = v1(N_o, basis)#; quiet=true)

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