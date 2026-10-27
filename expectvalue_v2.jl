# Version 3 of expectvalue_v2.jl with a --quiet (-q) flag. Everything else is unchanged.
#
# Usage:
#     julia -t <nthreads> expectvalue_v3.jl -f <state file or folder> -i <interaction file> [--decimal --n_orb <n_orb>] [-q]
#
# With --quiet, the only printout is the variational energy: one number per line, one line per state
# (for a folder, in the order g_0, g_1, ...). The state files are read with readwf(...; verbose=false), and the
# timing from @time is skipped. Error messages (missing files etc.) are still printed.
# The energy file is written in both modes.

include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("SphereED.jl")
include("SphereEDThreeBody.jl")
using .FQH_states
using .TwoBodyGet
using .ThreeBodyGet
using .ConstructManybodyMatrix
using LinearAlgebra
using ArgMacros
using Printf

# Variational energy of n-body Haldane pseudopotentials, computed from "removal amplitudes".
#
# For an n-tuple T = (t1<...<tn) of orbitals and a state ψ, let x_R(T) = <R| c_tn ... c_t1 |ψ>,
# where R runs over (N-n)-particle states. An n-body operator Σ_{T,T'} H[T,T'] c†_t1...c†_tn c_t'n...c_t'1
# then has <ψ|H|ψ> = Σ_R x_R† H x_R. H only connects tuples with the same index sum ("sector"), and in each
# sector it is written as H = Σ_k w_k W_k W_k†, so that <ψ|H|ψ> = Σ_k w_k Σ_R |W_k† x_R|².
# Each basis state of ψ is visited once per occupied n-tuple, so the cost grows as dim*N^n instead of dim².
#
# Orbital index i = 1..n_orb is bit (i-1) of the decimal basis, as in SphereED.jl.

struct Sector
    nb::Int                  # 2 or 3 body
    tidx::Matrix{Int}        # tidx[a,b] = row of W for the tuple starting with (a,b) (the last index is fixed by the sector)
    W::Matrix{ComplexF64}    # columns W_k
    w::Vector{Float64}       # weights w_k
end

# Two-body: W_v[(a,b)] = <(2s-m),M|a,b> from pp_matrix, w_v = c_v.
# This gives the same normalization as expectvalue.jl and sphere_FQHE_parallel.jl.
function twobody_sectors(n_orb::Int, v_list::Vector{Int32}, c_list::Vector{Float64})
    s = (n_orb-1)/2
    vmats = [pp_matrix(s, v) for v in v_list]
    sectors = Dict{Int,Sector}()
    for S in 3:(2n_orb-1)
        pairs = [(a,S-a) for a in 1:n_orb if a < S-a <= n_orb]
        tidx = zeros(Int, n_orb, n_orb)
        for (p,(a,b)) in enumerate(pairs)
            tidx[a,b] = p
        end
        W = zeros(ComplexF64, length(pairs), 0)
        w = Float64[]
        for (vmat, c) in zip(vmats, c_list)
            col = [vmat[a,b] for (a,b) in pairs]
            if any(!iszero, col)
                W = hcat(W, col)
                push!(w, c)
            end
        end
        if length(w) > 0
            sectors[S] = Sector(2, tidx, W, w)
        end
    end
    return sectors
end

# Three-body: H[T,T'] from get_threebody, diagonalized sector by sector.
# This gives the same normalization as sphere_FQHE_parallel_3body.jl (Vₘ with coefficient 1 is a projector).
function threebody_sectors(n_orb::Int, v_list::Vector{Int}, c_list::Vector{Float64})
    threeBody = get_threebody(n_orb, v_list, c_list)
    triples = ThreeBodyGet.triple_sectors(n_orb)
    H = Dict(S => zeros(ComplexF64, length(t), length(t)) for (S,t) in triples)
    index = Dict(t => k for (_,ts) in triples for (k,t) in enumerate(ts))
    for (h,(i,j,k,l,m,n)) in threeBody
        H[i+j+k][index[(i,j,k)], index[(l,m,n)]] += h
    end
    sectors = Dict{Int,Sector}()
    for (S,ts) in triples
        λ, U = eigen(Hermitian(H[S]))
        keep = findall(x -> abs(x) > 1e-12, λ)
        if isempty(keep) continue end
        tidx = zeros(Int, n_orb, n_orb)
        for (p,(a,b,_)) in enumerate(ts)
            tidx[a,b] = p
        end
        sectors[S] = Sector(3, tidx, U[:,keep], λ[keep])
    end
    return sectors
end

# Σ_k w_k Σ_R |W_k† x_R|² for one sector, accumulating W_k† x_R in a table indexed by R
function sector_energy(sec::Sector, S::Int, basis::Vector{Int64}, coefs::Vector{ComplexF64}, occ::Matrix{Int8})
    r = length(sec.w)
    Wc = conj.(sec.W)
    rowof = Dict{Int64,Int}()
    Y = ComplexF64[]
    N = size(occ,1)
    for (IInd, I) in enumerate(basis)
        ψ = coefs[IInd]
        if iszero(ψ) continue end
        for p in 1:N
            a = Int(occ[p,IInd])
            for q in (p+1):N
                b = Int(occ[q,IInd])
                if sec.nb == 2
                    if a+b != S continue end
                    J, sign = C(a, I, 1)
                    J, sign = C(b, J, sign)
                else
                    c = S-a-b
                    if c <= b break end # occ is sorted, so c only decreases as b increases
                    if (I>>(c-1))&1 == 0 continue end
                    J, sign = C(a, I, 1)
                    J, sign = C(b, J, sign)
                    J, sign = C(c, J, sign)
                end
                t = sec.tidx[a,b]
                n = length(rowof)
                row = get!(rowof, J, n+1)
                if row == n+1
                    append!(Y, zeros(ComplexF64, r))
                end
                amp = sign*ψ
                for k in 1:r
                    Y[(row-1)*r+k] += Wc[t,k]*amp
                end
            end
        end
    end
    E = 0.0
    for row in 1:length(rowof), k in 1:r
        E += sec.w[k]*abs2(Y[(row-1)*r+k])
    end
    return E
end

function nbody_energy(sectors::Dict{Int,Sector}, basis::Vector{Int64}, coefs::Vector{ComplexF64}, occ::Matrix{Int8})
    keys_list = collect(keys(sectors))
    if isempty(keys_list) || size(occ,1) < sectors[keys_list[1]].nb return 0.0 end
    energies = zeros(length(keys_list))
    Threads.@threads :dynamic for k in eachindex(keys_list)
        S = keys_list[k]
        energies[k] = sector_energy(sectors[S], S, basis, coefs, occ)
    end
    return sum(energies)
end

# Occupied orbitals of each basis state (sorted), one column per state
function occupation_table(basis::Vector{Int64}, n_orb::Int)
    N = count_ones(basis[1])
    occ = zeros(Int8, N, length(basis))
    for (IInd, I) in enumerate(basis)
        if count_ones(I) != N
            error("All basis states must have the same number of electrons.")
        end
        p = 0
        for i in 1:n_orb
            if (I>>(i-1))&1 == 1
                p += 1
                occ[p,IInd] = i
            end
        end
    end
    return occ
end

# Interaction file: "m c" for a two-body Vₘ, "m c 3" for a three-body Vₘ. Blank lines and lines starting with # are skipped.
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
    @inlinearguments begin
        @helpusage "julia -t <nthreads> expectvalue_v3.jl -f <state file or folder> -i <interaction file> [--decimal --n_orb <n_orb>] [-q]"
        @helpdescription """
        Variational energy of a given state with two- and three-body Haldane pseudopotentials on the sphere.
        Each line of the interaction file is 'm c' (two-body Vₘ with coefficient c) or 'm c 3' (three-body Vₘ).
        Two-body energies follow expectvalue.jl/sphere_FQHE_parallel.jl and three-body energies follow
        sphere_FQHE_parallel_3body.jl. The state is normalized before computing the energy."""
        @argumentrequired String fname "-f" "--file-name"
        @arghelp "State file, or a folder: then every g_* file in it is evaluated."
        @argumentrequired String intname "-i" "--interaction-file"
        @arghelp "Interaction file."
        @argumentflag decimal "--decimal"
        @arghelp "The state file(s) use the decimal basis format (requires --n_orb)."
        @argumentoptional Int N_o "--n_orb" "-o"
        @arghelp "Number of orbitals (only needed with --decimal)."
        @argumentflag quiet "--quiet" "-q"
        @arghelp "Only print the variational energy."
    end
    say(x...) = quiet || println(x...) # println that is silent with --quiet
    say("Calculate two- and three-body pseudopotential variational energy of a given state.")

    if decimal && N_o == nothing
        println("For decimal format, the number of orbital must be specified with '--n_orb'.")
        return
    end

    if isdir(fname)
        files = filter(f -> occursin(r"^g_\d+$", f), readdir(fname))
        sort!(files, by = f -> parse(Int, f[3:end]))
        files = joinpath.(fname, files)
        if isempty(files)
            println("No g_* files found in '$(fname)'. Terminating.")
            return false
        end
    elseif isfile(fname)
        files = [fname]
    else
        println("State file '$(fname)' not found. Terminating.")
        return false
    end

    if !isfile(intname)
        print("Interaction file '$(intname)' not found. Terminating.")
        return false
    end
    v2, c2, v3, c3 = read_interaction(intname)
    say("Two-body:   m = $(v2), c = $(c2)")
    say("Three-body: m = $(v3), c = $(c3)")

    sectors2 = nothing
    sectors3 = nothing
    n_orb_prev = 0
    results = Tuple{String,Float64}[]
    for file in files
        say("--------")
        say("State: $(file)")
        state = decimal ? readwfdec(file, N_o; verbose=!quiet) : readwf(file; verbose=!quiet)
        n_orb = length(state.basis[1])
        basis = decimalbasis(state)
        coefs = ComplexF64.(state.coef)
        say("$(n_orb) orbitals, $(count(state.basis[1])) electrons.")

        # Matrix elements only depend on n_orb, so they are reused for every state in a folder
        if n_orb != n_orb_prev
            sectors2 = twobody_sectors(n_orb, v2, c2)
            sectors3 = threebody_sectors(n_orb, v3, c3)
            n_orb_prev = n_orb
        end

        norm2 = sum(abs2, coefs)
        say("<ψ|ψ> = $(norm2)")
        occ = occupation_table(basis, n_orb)
        energies() = (isempty(v2) ? 0.0 : nbody_energy(sectors2, basis, coefs, occ)/norm2,
                      isempty(v3) ? 0.0 : nbody_energy(sectors3, basis, coefs, occ)/norm2)
        E2, E3 = quiet ? energies() : @time energies()
        if !isempty(v2) && !isempty(v3)
            say("Two-body energy   = $(E2)")
            say("Three-body energy = $(E3)")
        end
        quiet ? println(E2+E3) : println("The energy is $(E2+E3)")
        push!(results, (file, E2+E3))
    end

    say("--------")
    intbase = splitext(basename(intname))[1]
    outname = isdir(fname) ? joinpath(fname, "energy_$(intbase).txt") : "$(fname)_$(intbase)_energy.txt"
    open(outname, "w") do f
        for (file, E) in results
            if isdir(fname)
                write(f, "$(basename(file)) $(E)\n")
            else
                write(f, "$(E)\n")
            end
        end
    end
    say("Saved energy to $(outname)")
end

# @time is skipped with --quiet (the flag is parsed inside main, so ARGS is checked here)
if any(a -> a in ("-q", "--quiet"), ARGS)
    main()
else
    @time main()
end
