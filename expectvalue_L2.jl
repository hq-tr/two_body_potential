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
    println("Diagonalize L⁺L⁻ matrix.")
    println("INSTRUCTION: Place all the basis states in a subdirectory. They will be read automatically.")
    println("Please ensure that subdirectory contains only the basis states and nothing else.")
    println()
    println("IMPORTANT: It is assumed that all the basis states have the same monomials listed in the same order. The code will not work if that is not the case.")
    println()
    println("Input sub-director file name: ")

    dirname = readline()
    if !isdir(dirname)
        println("$(dirname) is not a directory! EXITING.")
        return
    end
    f_opt = 0

    while !(f_opt in [1,2])
        println("Are your vectors in binary (1) or decimal (2) format?")
        f_opt = parse(Int,readline())
        if !(f_opt in [1,2])
            println("Please input only 1 or 2")
        end
    end

    if f_opt == 2
        println("How many orbital?")
        No = parse(Int,readline())
    end
    all_states = Vector{Float64}[]
    basis  = BitVector[]
    n_files = 0
    for filename in readdir(dirname)
        print("Reading $(filename):    ")
        if f_opt == 1
            readstate = readwf("$(dirname)/$(filename)")
        else
            readstate = readwfdec("$(dirname)/$(filename)",No)
        end
        if n_files == 0 
            for vec in readstate.basis
                push!(basis,vec)
            end
        end
        push!(all_states,readstate.coef)
        n_files += 1
    end

    println("$(n_files) file(s) found.")

    #basis = all_basis[0]
    dim   = length(basis)
    coefs = reduce(hcat,all_states) # this is a matrix each of whose columns is a vector that corresponds to a basis state

    N_o = length(basis[1])
    println("$(N_o) orbitals.")

    #@time basis, dim = getbasis(filewf, N_o, N_e)

    println("--------")
    println("Constructing the Hamiltonian")

    @time H_matrix = L⁺L⁻(basis)
    @time H_ED     = transpose(coefs)*H_matrix*coefs

    eval = eigvals(H_ED)
    evec = eigvecs(H_ED)

    if !isdir("$(dirname)_eigen") 
        mkdir("$(dirname)_eigen") 
    end

    open("$(dirname)_eigen/L2_eigenvalues.txt", "w+") do f
        for ev in eval 
            write(f,"$(ev)\n") 
        end
    end

    eigen_coefs = coefs * evec
    for i in 1:n_files
        eigen_state = FQH_state(basis, eigen_coefs[:,i])
        printwf(eigen_state;fname="$(dirname)_eigen/L2_eigen_$(i-1)")
    end
    #@time H_matrix = two_body(N_o, basis, v_list, c_list)

    #println("--------")
    #println(size(coefs))
    #println(size(H_matrix))

    #@time ϵ = coefs' * H_matrix * coefs
    println("Done.")


#    println(transpose(ϕ))

end

@time main()