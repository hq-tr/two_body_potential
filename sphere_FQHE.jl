include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("/home/trung/_qhe-julia/HilbertSpace.jl")
include("v1_sphere.jl")
using .FQH_states
using .PseudoPotential
using .HilbertSpaceGenerator
using LinearAlgebra
using Arpack
using SparseArrays
using ArgMacros
using BenchmarkTools
using Printf

function main()
    # ================================ READ USER INPUT ================================
    @inlinearguments begin
        @argumentdefault String "" fname "-f" "--filename"
        @argumentdefault Int 5 k "-n" "--nev"
        @argumentdefault String "" intname "-i" "--interaction-file"
        @argumentoptional Int n_el "-e" "--n_el"
        @argumentoptional Int n_orb "-o" "--n_orb"
        @argumentoptional Float64 L_z "-Lz" "--Lz"
    end

    println("============================================================")
    println("      FULL-ED OF TWO-BODY INTERACTION ON THE SPHERE")
    println("============================================================")

    if n_el != nothing && n_orb != nothing
        if L_z == nothing
            println("Generating a basis with $(n_el) electrons and $(n_orb) orbitals.")
            basis = fullhilbertspace(n_el,n_orb)
            outname = @sprintf "%ie_%io" n_el n_orb
        else
            println("Generating a basis with $(n_el) electrons and $(n_orb) orbitals in the Lz=$(L_z) sector.")
            basis = fullhilbertspace(n_el,n_orb,L_z)
            outname = @sprintf "%ie_%io_%.1f" n_el n_orb L_z
        end
    elseif length(fname) > 0
        println("Reading basis vectors from [$(fname)]")
        
        basis = readwf(fname).basis

        outname = fname
    else
        println()
        println("WARNING: No input or incomplete input was specified. The program will now terminating.")
        println("Run the program with '-h' or '--help' tag to view possible arguments.")
        println()
        return
    end

    N_o = length(basis[1])
    d   = length(basis) # Dimension

    v_list = Int32[]
    c_list = Float64[]

    if length(intname) == 0
        println("Input m for Vₘ and the corresponding coefficient. ")
        println("Each pp term takes one line, with two numbers separated by a space.")
        println("Put a 0 to end")
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
    else
        println("Reading interaction from $(intname).")
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
    end

    #@time basis, dim = getbasis(filewf, N_o, N_e)

    # ======================== CONSTRUCT AND DIAGONALIZE HAMILTONIAN ======================
    println("--------")
    println("Constructing the Hamiltonian")

    @time H_matrix = two_body(N_o, basis, v_list, c_list)#; quiet=true)

    display(H_matrix)

    println("--------")

    println("Diagonalizing with ARPACK")

    @time λ, ϕ = eigs(H_matrix, nev=k,which=:SM)

    #display(ϕ)

    # ====================== SAVE GROUND STATE =======================
    println("Eigenvalues = ")
    display(real.(λ))

    println("--------")
    dirname = "$(outname)_$(intname)_out"

    if !isdir(dirname) mkdir(dirname) end

    open("$(dirname)/eigen.txt","w+") do f
        for value in λ
            write(f,"$(value)\n")
        end
    end

    for i in 1:k
        gs_coef = ϕ[:,i]
        #println(length(gs_coef))
        ground_state = FQH_state(basis, gs_coef)
        printwf(ground_state;fname="$(dirname)/g_$(fname)_$(i-1)")
    end


    #println("Saved ground state as g_$(fname)_0.")


#    println(transpose(ϕ))

end

@time main()