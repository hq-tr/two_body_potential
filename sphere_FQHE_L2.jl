include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("v1_sphere.jl")
using .FQH_states
using .PseudoPotential
using LinearAlgebra
using Arpack
using SparseArrays
using ArgMacros
using BenchmarkTools

function main()
    # ================================ READ USER INPUT ================================
    @inlinearguments begin
        @argumentrequired String fname "-f" "--filename"
        @argumentrequired Int k "-n" "--nev"
    end
    
    basis = readwf(fname).basis

    N_o = length(basis[1])
    d   = length(basis) # Dimension


    println("$(N_o) orbitals.")

    # ======================== CONSTRUCT AND DIAGONALIZE HAMILTONIAN ======================
    println("--------")
    println("Constructing the Hamiltonian")

    @time H_matrix = L⁺L⁻(basis)#; quiet=true)

    display(H_matrix)

    println("--------")

    println("Diagonalizing with ARPACK")

    @time λ, ϕ = eigs(H_matrix, nev=k,which=:SM)


    println("Eigenvalues = ")
    if k < 15
        display(real.(λ))
    else
        display(real.(λ[1:15]))
        println("⋮")
    end

    # ====================== SAVE ZERO-EIGENVALUE EIGENSTATES =======================
    println("--------")
    i = 1
    tol = 1e-10 # tolerance

    dirname = fname * "_groundstates" # All ground states will be saved in this folder
    if !isdir(dirname)
        mkdir(dirname)
    end

    # Save the eigenvalues
    open("$(dirname)/eigen.txt", "w+") do f
        for ε in λ write(f,"$ε\n") end
    end

    while λ[i] < tol 
        gs_coef = ϕ[:,1] # these are the co-efficients
        gs = FQH_state(basis,gs_coef)
        printwf(gs;fname="$(dirname)/g_$(i-1)")
        i += 1
    end


#    println(transpose(ϕ))

end

@time main()