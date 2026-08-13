include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("/home/trung/_qhe-julia/HilbertSpace.jl")
include("/home/trung/_qhe-julia/Potentials.jl")
include("/home/trung/_qhe-julia/Misc.jl")
include("SphereED.jl")
include("PseudoPotentials.jl")
using .FQH_states
using .PseudoPotential
using .HilbertSpaceGenerator
using .Potentials
using .MiscRoutine
using .OneBodyGet
using .TwoBodyGet
using .ConstructManybodyMatrix
using LinearAlgebra
using Arpack
using SparseArrays
using ArgMacros
using BenchmarkTools
using Printf
using Combinatorics

function get_rows_and_cols(basis::Vector{T} where T<:Integer,n_orb::Integer;onebody=true,twobody=true)
    if !onebody && !twobody return UInt64[],UInt64[] end

    d = length(basis)

    #for b in basis
    #   println(bitstring(b))
    #end

    # Initialize the rows and cols vectors with diagonal terms
    rows = UInt64.(1:d)
    cols = UInt64.(1:d)

    # For every element in basis, find all other elements with which it gives non-zero matrix elements
    for (index1,basis1) in enumerate(basis)
        #print("\r   Checking basis $(index1) of $d\t\t\t")
        basis1_dex = dec2dex(basis1)

        # One-body c†_m2 c_m1
        if onebody
            for m1 in basis1_dex
                for m2 in (m1+1):(n_orb-1) #only check m2 > m1
                    
                    if m2 in basis1_dex continue end # Make sure m2 is not an electron

                    basis2 = basis1 - 2^m1 + 2^m2
                    index2 = searchsortedfirst(basis,basis2)

                    push!(rows,index1)
                    push!(cols,index2)
                    push!(rows,index2)
                    push!(cols,index1)

                end
            end
        end

        # Two-body c†_m3 c†_m4 c_m1 c_m2 (m3 < m1 < m2 < m4)
        if twobody
            for (m1,m2) in combinations(basis1_dex,2)
                for m4 in (m2+1):min(n_orb,m2+m1)
                    if m4 in basis1_dex continue end
                    m3 = m1+m2-m4
                    if m3 in basis1_dex continue end
                    if m3 > n_orb continue end

                    basis2 = basis1 - 2^m1 - 2^m2 + 2^m3 + 2^m4
                    index2 = searchsortedfirst(basis,basis2)
                    if index2 > d continue end
                    if basis[index2] != basis2 continue end

                    push!(rows,index1)
                    push!(cols,index2)
                    push!(rows,index2)
                    push!(cols,index1)
                end
            end
        end
    end
    return rows,cols
end

function main()
    # ================================ READ USER INPUT ================================
    @inlinearguments begin
        @argumentoptional String fname "-f" "--filename"
        @argumentdefault Int 5 k "-n" "--nev"
        @argumentdefault String "" intname "-i" "--interaction-file"
        @argumentdefault String "" pinname "-p" "--pin-potential-file"
        @argumentrequired Int n_el "-e" "--n_el"
        @argumentrequired Int n_orb "-o" "--n_orb"
        @argumentoptional Float64 L_z "-Lz" "--Lz"
        @argumentdefault Float64 3.141592653589793 angle_multiplier "--angle-multiplier" 
        @argumentflag quiet "--quiet"
        @argumentdefault Int 1 numzones1bd "-z" "--num-zones-1bdy"
        @argumentdefault Int 1 numzones2bd "-Z" "--num-zones-2bdy"
    end

    println("============================================================")
    println("      FULL-ED OF TWO-BODY INTERACTION ON THE SPHERE")
    println("============================================================")

    # Reading basis input
    if fname != nothing
        println("Reading basis vectors from [$(fname)]")
        state = readwf(fname)
        N_o   = length(state.basis[1])
        N_e   = count(state.basis[1])
        @assert((n_el==N_e) && (n_orb==N_o),"Input number of electrons and number of orbitals don't match the file.")
        basis = decimalbasis(state;sortorder=:ascending) # The algorithm requires the basis in a sorted order
        outname = fname
    else
        if L_z == nothing
            println("Generating a basis with $(n_el) electrons and $(n_orb) orbitals (all Lz sectors).")
            basis = fullhilbertspace(n_el,n_orb;output_type="Decimal")
            outname = @sprintf "%ie_%io" n_el n_orb
        else
            println("Generating a basis with $(n_el) electrons and $(n_orb) orbitals in the Lz=$(L_z) sector.")
            basis = fullhilbertspace(n_el,n_orb,L_z;output_type="Decimal")
            outname = @sprintf "%ie_%io_%.1f" n_el n_orb L_z
        end
    end

    N_o = n_orb
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
    end

    # Reading pin input
    θ_list = Float64[]
    ϕ_list = Float64[]
    V_list = Float64[]
    if length(pinname) > 0
        if isfile(pinname)
            println("Reading pinning potential data from $(pinname).")
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
    println("$npins pin(s).")
    if npins > 0
        pinappendname = "$(pinname)_"
    else
        pinappendname = ""
    end
    #@time basis, dim = getbasis(filewf, N_o, N_e)

    # ======================== CONSTRUCT AND DIAGONALIZE HAMILTONIAN ======================
    println("--------")
    println("Constructing the Hamiltonian")   


    # Pre-allocate memory for the matrix
    if npins == 0 # only two-body
        rows, cols = get_rows_and_cols(basis,n_orb;onebody=false)
    elseif length(v_list) == 0 # only one-body
        rows, cols = get_rows_and_cols(basis,n_orb;twobody=false)
    else # both one-body and two-body
        rows, cols = get_rows_and_cols(basis,n_orb)
    end

    nnz = length(rows)

    H_matrix   = spzeros(ComplexF64,rows,cols,d,d) # requires Julia 1.10 or newer
        
    rows = nothing
    cols = nothing
    GC.gc() # free up memory

    println("Pre-allocated a sparse matrix with $(nnz) non-zero elements")

    # Construct two-body
    if length(v_list) > 0
        println("\nConstructing the two-body matrix")
        twoBody = get_twobody(n_orb, v_list::Vector{Int32}, c_list::Vector{Float64})
        @time begin
            for subzone = 1:numzones2bd
                if !quiet
                    print("\r   Working on subzone $(subzone)/$(numzones2bd)\t\t\t")
                end
                calcV!(H_matrix,twoBody,basis,subzone;num_zones=numzones2bd);
            end
            println("Done!")
        end 
    end
    

    # Construct one-body
    if npins > 0
        println("\nConstructing the one-body matrix")
        pinsize_list = ones(Int,npins)
        oneBody = get_oneBody_widebump(n_orb,θ_list, ϕ_list, pinsize_list)
        @time begin
            for subzone = 1:numzones1bd
                print("\r   Working on subzone $(subzone)/$(numzones1bd)\t\t\t")
                calcT!(H_matrix,oneBody,basis,subzone,1.0;num_zones=numzones1bd)
            end
            println("Done!")
        end
    end

    # Display matrix
    if !quiet
        println("\nH matrix = ")
        if d < 20
            display(Matrix(H_matrix))
        elseif d< 1000
            display(H_matrix) # display sparse pattern
            # This is quite memory-intensive, so it is only done for smaller matrices
        else
            println(summary(H_matrix)) 
            println()
        end
    end

    println("--------")
    println("Diagonalizing with ARPACK")

    @time λ, ϕ = eigs(H_matrix, nev=k,which=:SM)

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
        ground_state = FQH_state(basis, gs_coef, n_orb)
        printwf(ground_state;fname="$(dirname)/g_$(i-1)")
    end


    #println("Saved ground state as g_$(fname)_0.")


#    println(transpose(ϕ))

end

@time main()