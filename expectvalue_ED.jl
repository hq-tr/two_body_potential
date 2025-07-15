include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("/home/trung/two_body_potential/v1_sphere.jl")
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
        #@argumentrequired Int dim "-d" "--dimension"
        @argumentrequired String intname "-i" "--interaction-file"
    end

println("============")
    println("IMPORTANT: all the input states must contain the same monomial basis listed in the same order.")
    println("The input files must be named in the following format: <name root><index>")
    println("where <index> is an integer from 0 to dim-1.")
    println("------")
    println("Reading files with name root '$(fname)'")

    if isfile(fname * "0")
        state = readwf(fname * "0")
        basis = state.basis
        h_dim = length(basis)
        coefs = Float64[]
    else
        println("File name '$(fname)0' not found. Terminating.")
        return false
    end
    d=0
    while (isfile(fname * string(d)))
        state = readwf(fname * string(d))
        append!(coefs,identity.(state.coef))
        d+=1
    end

    println("Found $d file(s).")
    println("-----------------------")

    coefs = reshape(coefs,(h_dim,d))

    N_o = length(basis[1])
    println("$(N_o) orbitals.")


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

    if d == 1
        @time ϵ = coefs' * H_matrix * coefs
        println("The energy is $ϵ")
        println("--------")
    else
        @time mat = coefs' * H_matrix * coefs
        ϵ = eigvals(mat)
        println("The eigenvalues are")
        display(ϵ)

        # ======= Save the outputs
        # Create directory
        dirname = "$(fname)_$(intname)"
        if !isdir(dirname) mkdir(dirname) end
        
        # Save the eigenvalues
        open("$(dirname)/eigen.txt","w+") do f
            for value in ϵ
                write(f,"$(real(value))\n")
            end
        end

        # Save the eigenstates
        vecs = eigvecs(mat)

        for i in 1:d
            state_coef = coefs * vecs[:,i]
            state = FQH_state(basis, state_coef)
            printwf(state;fname="$(dirname)/g_$(i-1)")
            i+=1
        end
    end
    



    return true
#    println(transpose(ϕ))

end

@time main()