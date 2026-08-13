include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("/home/trung/_qhe-julia/HilbertSpace.jl")
include("/home/trung/_qhe-julia/Potentials.jl")
include("PseudoPotentials.jl")
using .FQH_states
using .PseudoPotential
using .HilbertSpaceGenerator
using .Potentials
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
        @argumentdefault String "" pinname "-p" "--pin-potential-file"
        @argumentoptional Int n_el "-e" "--n_el"
        @argumentoptional Int n_orb "-o" "--n_orb"
        @argumentoptional Float64 L_z "-Lz" "--Lz"
        @argumentdefault Float64 3.141592653589793 angle_multiplier "--angle-multiplier" 
        @argumentflag verbose "--verbose"
    end

    if verbose
        println("============================================================")
        println("      FULL-ED OF TWO-BODY INTERACTION ON THE SPHERE")
        println("============================================================")
    end
    # Reading basis input
    if n_el != nothing && n_orb != nothing
        if L_z == nothing
            println("Generating a basis with $(n_el) electrons and $(n_orb) orbitals (all Lz sectors).")
            basis = fullhilbertspace(n_el,n_orb)
            outname = @sprintf "%ie_%io" n_el n_orb
        else
            println("Generating a basis with $(n_el) electrons and $(n_orb) orbitals in the Lz=$(L_z) sector.")
            basis = fullhilbertspace(n_el,n_orb,L_z)
            outname = @sprintf "%ie_%io_%.1f" n_el n_orb L_z
        end
    elseif length(fname) > 0
        println("Reading basis vectors from [$(fname)]")
        
        basis = readwf(fname;verbose=verbose).basis

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

    # Reading two-body interaction input
    v_list = Int32[]
    c_list = Float64[]


    if intname != "none"
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
            if verbose
                println("Reading interaction from $(intname).")
            end
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
    end

    # Reading pin input
    θ_list = Float64[]
    ϕ_list = Float64[]
    V_list = Float64[]
    if length(pinname) > 0
        if isfile(pinname)
            if verbose
                println("Reading pinning potential data from $(pinname).")
            end
            open(pinname) do f
                for line in map(s->split(s),readlines(f))
                    append!(θ_list,parse(Float64,line[1])*angle_multiplier)
                    append!(ϕ_list,parse(Float64,line[2])*angle_multiplier)
                    append!(V_list,parse(Float64,line[3]))
                end
            end
        else
            println("Pinning potential data file '$(pinname)' not found. Terminating.")
        end
    end
    npins = length(θ_list)
    if verbose
        println("$npins pin(s).")
    end
    if npins > 0
        pinappendname = "$(pinname)_"
    else
        pinappendname = ""
    end
    #@time basis, dim = getbasis(filewf, N_o, N_e)

    # ======================== CONSTRUCT AND DIAGONALIZE HAMILTONIAN ======================
    if verbose
        println("--------")
        println("Constructing the Hamiltonian")
    end

    if length(v_list) > 0
        @time H_matrix = two_body_sphere(N_o, basis, v_list, c_list)#; quiet=true)
    else
        @time H_matrix = spzeros(ComplexF64,(d,d))
    end

    @time for i in 1:npins
        H_matrix += sphere_point_matrix(basis, θ_list[i], ϕ_list[i], V_list[i])
    end

    if verbose
        println("Hamiltonian matrix = ")
        display(H_matrix) 
    end
    open("matrix.dat","w+") do f
        matnz = findnz(H_matrix)
        for i in 1:length(matnz[1])
            write(f,"$(matnz[1][i]),$(matnz[2][i]),$(matnz[3][i])\n")
        end
    end

    if verbose
        println("--------")
        println("Diagonalizing with ARPACK")
    end

    @time λ, ϕ = eigs(H_matrix, nev=k,which=:SR)

    #display(ϕ)

    # ====================== SAVE GROUND STATE =======================
    println("Eigenvalues = ")
    display(real.(λ))

    println("--------")
    dirname = "$(outname)_$(intname)_$(pinappendname)out"

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
        printwf(ground_state;fname="$(dirname)/g_$(i-1)")
    end


    #println("Saved ground state as g_$(fname)_0.")


#    println(transpose(ϕ))

end

@time main()