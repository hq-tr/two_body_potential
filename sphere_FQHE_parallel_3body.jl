# Version 2 of sphere_FQHE_parallel_3body.jl with a single interaction file, in the same format as expectvalue_v2.jl.
#
# Usage:
#     julia -t <nthreads> sphere_FQHE_parallel_3body_v2.jl -e <n_el> -o <n_orb> -i <interaction file> [options]
#
# The -t (--threebody-file) option is removed. Both two- and three-body pseudopotentials are read from the
# file given with -i (now required). Each line of the interaction file is
#     m c      two-body Vₘ with coefficient c
#     m c 3    three-body Vₘ with coefficient c
# Blank lines and lines starting with # are skipped.
# Only the input handling changes; the Hamiltonian, diagonalization and normalization are as in
# sphere_FQHE_parallel_3body.jl.
#
# The output directory is <basis name>_<interaction file>_[<pin file>_]out.

include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("/home/trung/_qhe-julia/HilbertSpace.jl")
include("/home/trung/_qhe-julia/Potentials.jl")
include("/home/trung/_qhe-julia/Misc.jl")
include("SphereED.jl")
include("SphereEDThreeBody.jl")
include("PseudoPotentials.jl")
using .FQH_states
using .PseudoPotential
using .HilbertSpaceGenerator
using .Potentials
using .MiscRoutine
using .OneBodyGet
using .TwoBodyGet
using .ThreeBodyGet
using .ConstructManybodyMatrix
using LinearAlgebra
using Arpack
using SparseArrays
using ArgMacros
using Printf
using Combinatorics

# Sparsity pattern of the Hamiltonian.
# Three-body: c†_n1 c†_n2 c†_n3 c_m1 c_m2 c_m3 with n1+n2+n3 = m1+m2+m3. This includes every two-body hop
# (with the third particle as a spectator), so two-body terms are covered as long as n_el >= 3.
# One-body: c†_m2 c_m1 for any pair, needed only when there are pins.
# Orbitals here are 0-based (bit m of a basis state is orbital m).
function get_rows_and_cols_3body(basis::Vector{T} where T<:Integer,n_orb::Integer;onebody=true,threebody=true)
    if !onebody && !threebody return UInt64[],UInt64[] end

    d = length(basis)
    connected = Vector{Vector{UInt64}}(undef,d)

    Threads.@threads for index1 in 1:d
        basis1 = basis[index1]
        basis1_dex = [m for m in 0:(n_orb-1) if (basis1>>m)&1 == 1]
        found = UInt64[index1] # diagonal term

        # One-body c†_m2 c_m1
        if onebody
            for m1 in basis1_dex
                for m2 in 0:(n_orb-1)
                    if (basis1>>m2)&1 == 1 continue end # Make sure m2 is not an electron
                    basis2 = basis1 - (1<<m1) + (1<<m2)
                    index2 = searchsortedfirst(basis,basis2)
                    if index2 > d || basis[index2] != basis2 continue end
                    push!(found,index2)
                end
            end
        end

        # Three-body: remove (m1,m2,m3), then put back (n1<n2<n3) with the same total
        if threebody
            for (m1,m2,m3) in combinations(basis1_dex,3)
                rest = basis1 - (1<<m1) - (1<<m2) - (1<<m3)
                msum = m1+m2+m3
                for n1 in 0:(n_orb-1)
                    if (rest>>n1)&1 == 1 continue end
                    for n2 in (n1+1):(n_orb-1)
                        n3 = msum-n1-n2
                        if n3 <= n2 break end
                        if n3 > n_orb-1 continue end
                        if (rest>>n2)&1 == 1 || (rest>>n3)&1 == 1 continue end
                        basis2 = rest + (1<<n1) + (1<<n2) + (1<<n3)
                        index2 = searchsortedfirst(basis,basis2)
                        if index2 > d || basis[index2] != basis2 continue end
                        push!(found,index2)
                    end
                end
            end
        end
        connected[index1] = unique!(sort!(found))
    end

    # The pattern is symmetric, so collecting every connection of every row gives both (i,j) and (j,i)
    rows = Vector{UInt64}(undef,sum(length,connected))
    cols = similar(rows)
    pos = 0
    for (index1,c) in enumerate(connected)
        for index2 in c
            pos += 1
            rows[pos] = index2
            cols[pos] = index1
        end
    end
    return rows,cols
end

# Interaction file: "m c" for a two-body Vₘ, "m c 3" for a three-body Vₘ. Blank lines and lines starting with # are skipped.
# (same as read_interaction in expectvalue_v2.jl)
function read_interaction(intname::String)
    v2 = Int32[]; c2 = Float64[]
    v3 = Int[];   c3 = Float64[]
    for (n,line) in enumerate(readlines(intname))
        data = split(line)
        if isempty(data) || startswith(data[1], "#") continue end
        if length(data) < 2
            error("Line $n of $(intname) should be 'm coefficient [n_body]'.")
        end
        nb = length(data) >= 3 ? parse(Int, data[3]) : 2
        if nb == 2
            push!(v2, parse(Int32, data[1])); push!(c2, parse(Float64, data[2]))
        elseif nb == 3
            push!(v3, parse(Int, data[1])); push!(c3, parse(Float64, data[2]))
        else
            error("Line $n of $(intname): only 2- and 3-body interactions are supported.")
        end
    end
    return v2, c2, v3, c3
end

function main()
    # ================================ READ USER INPUT ================================
    @inlinearguments begin
        @helpusage "julia -t <nthreads> sphere_FQHE_parallel_3body_v2.jl -e <n_el> -o <n_orb> -i <interaction file> [options]"
        @helpdescription """
        Full ED of three-body Haldane pseudopotentials on the sphere (fermions), optionally together
        with two-body pseudopotentials and one-body potential pins.
        Each line of the interaction file is 'm c' (two-body Vₘ with coefficient c) or 'm c 3' (three-body Vₘ).
        Eigenvalues are in the same units as sphere_FQHE_parallel.jl, with each three-body Vₘ the
        projector onto three-particle states of total angular momentum 3S-m."""
        @argumentoptional String fname "-f" "--filename"
        @arghelp "Vector file whose basis is used. Overrides --Lz."
        @argumentdefault Int 5 k "-n" "--nev"
        @arghelp "Number of eigenvalues/eigenvectors to compute (default 5)."
        @argumentrequired String intname "-i" "--interaction-file"
        @arghelp "Interaction file with lines 'm c' (two-body) and 'm c 3' (three-body)."
        @argumentdefault String "" pinname "-p" "--pin-potential-file"
        @arghelp "Pin data file: polar angle, azimuthal angle, strength on each line."
        @argumentrequired Int n_el "-e" "--n_el"
        @arghelp "Number of electrons."
        @argumentrequired Int n_orb "-o" "--n_orb"
        @arghelp "Number of orbitals."
        @argumentoptional Float64 L_z "-Lz" "--Lz"
        @arghelp "Lz sector. All sectors are used if omitted."
        @argumentdefault Float64 3.141592653589793 angle_multiplier "--angle-multiplier"
        @arghelp "Multiplier for the pin angles (default π)."
        @argumentflag quiet "--quiet"
        @arghelp "Do not print the Hamiltonian or progress."
        @argumentdefault Int 1 numzones1bd "-z" "--num-zones-1bdy"
        @arghelp "Number of subzones for the one-body part (more = less memory, slower)."
        @argumentdefault Int 1 numzones2bd "-Z" "--num-zones-2bdy"
        @arghelp "Number of subzones for the two-body part."
        @argumentdefault Int 1 numzones3bd "-Y" "--num-zones-3bdy"
        @arghelp "Number of subzones for the three-body part."
    end

    println("============================================================")
    println("      FULL-ED OF THREE-BODY INTERACTION ON THE SPHERE")
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

    d   = length(basis) # Dimension
    println("Dimension = $d")
    if d == 0
        println("Empty Hilbert space. Terminating.")
        return false
    end

    # Reading interaction input (two- and three-body)
    if !isfile(intname)
        print("Interaction file '$(intname)' not found. Terminating.")
        return false
    end
    println("Reading interaction from $(intname).")
    v_list, c_list, v3_list, c3_list = read_interaction(intname)
    println("Two-body:   m = $(v_list), c = $(c_list)")
    println("Three-body: m = $(v3_list), c = $(c3_list)")

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
            return false
        end
    end
    npins = length(θ_list)
    println("$npins pin(s).")
    if npins > 0
        pinappendname = "$(pinname)_"
    else
        pinappendname = ""
    end
    intappendname = "$(intname)_"

    # ======================== CONSTRUCT AND DIAGONALIZE HAMILTONIAN ======================
    println("--------")
    println("Constructing the Hamiltonian")

    # Pre-allocate memory for the matrix. The three-body pattern also covers the two-body one when n_el >= 3.
    has_manybody = length(v3_list) > 0 || length(v_list) > 0
    if length(v_list) > 0 && n_el < 3
        rows, cols = get_rows_and_cols_3body(basis,n_orb;onebody=true,threebody=false)
        println("Fewer than 3 electrons: the pre-allocated pattern does not cover two-body terms.")
    else
        @time rows, cols = get_rows_and_cols_3body(basis,n_orb;onebody=(npins>0),threebody=has_manybody)
    end

    nnz = length(rows)

    H_matrix   = spzeros(ComplexF64,rows,cols,d,d) # requires Julia 1.10 or newer

    rows = nothing
    cols = nothing
    GC.gc() # free up memory

    println("Pre-allocated a sparse matrix with $(nnz) non-zero elements")

    # Construct three-body
    # The two-body list from get_twobody counts each term twice and the eigenvalues are halved at the end
    # (as in sphere_FQHE_parallel.jl), so the three-body coefficients are doubled here to compensate.
    if length(v3_list) > 0
        println("\nConstructing the three-body matrix")
        threeBody = get_threebody(n_orb, v3_list, 2 .* c3_list)
        println("   $(length(threeBody)) three-body matrix elements")
        @time begin
            for subzone = 1:numzones3bd
                if !quiet
                    print("\r   Working on subzone $(subzone)/$(numzones3bd)\t\t\t")
                end
                calcV3!(H_matrix,threeBody,basis,subzone;num_zones=numzones3bd);
            end
            println("Done!")
        end
        threeBody = nothing
        GC.gc()
    end

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
                if !quiet
                    print("\r   Working on subzone $(subzone)/$(numzones1bd)\t\t\t")
                end
                calcT!(H_matrix,oneBody,basis,subzone,1.0;num_zones=numzones1bd)
            end
            println("Done!")
        end
    end

    if nnz != SparseArrays.nnz(H_matrix)
        println("Warning: $(SparseArrays.nnz(H_matrix)-nnz) matrix elements were outside the pre-allocated pattern.")
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

    k = min(k, d-1)
    @time λ, ϕ = eigs(H_matrix, nev=k,which=:SM)
    λ = λ ./ 2 # to be consistent with the single-core version

    # ====================== SAVE EIGENSTATES =======================
    println("Eigenvalues = ")
    display(real.(λ))

    println("--------")
    dirname = "$(outname)_$(intappendname)$(pinappendname)out"

    if !isdir(dirname) mkdir(dirname) end

    open("$(dirname)/eigen.txt","w+") do f
        for value in λ
            write(f,"$(real(value))\n")
        end
    end

    for i in 1:k
        gs_coef = ϕ[:,i]
        ground_state = FQH_state(basis, gs_coef, n_orb)
        printwf(ground_state;fname="$(dirname)/g_$(i-1)")
    end
    println("Saved eigenstates to $(dirname)/")
end

@time main()
