module ThreeBodyGet
using LinearAlgebra
# Three-body Haldane pseudopotentials on the sphere (fermions).
#
# Conventions follow TwoBodyGet: orbital index i = 1..n_orb corresponds to Lz = s-i+1, s = (n_orb-1)/2.
# The three-body pseudopotential V_m^(3) is the projector onto three-particle states with total angular
# momentum L = 3s-m (summed over all multiplets if the L subspace is degenerate, e.g. m >= 9).
#
# Each output element (h,(i,j,k,l,m,n)) has i<j<k and l<m<n and stands for h*|i,j,k><l,m,n|, where
# |i,j,k> = c†_i c†_j c†_k|0>. In second quantization this is the term
#     h * c†_i c†_j c†_k c_n c_m c_l
# Summing these terms over the list gives the full Hamiltonian (no extra factor, no double counting).

# All ordered triples (a<b<c) of orbital indices, grouped by index sum a+b+c (i.e. by total Lz)
function triple_sectors(n_orb::Int)
    sectors = Dict{Int, Vector{NTuple{3,Int}}}()
    for a = 1:n_orb, b = a+1:n_orb, c = b+1:n_orb
        push!(get!(sectors, a+b+c, NTuple{3,Int}[]), (a,b,c))
    end
    return sectors
end

# Matrix of L- = Σ_i sqrt(s(s+1)-m_i(m_i-1)) c†_{i+1} c_i from sector `from` to sector `to`
function lowering_matrix(s::Float64, from::Vector{NTuple{3,Int}}, to::Vector{NTuple{3,Int}})
    n_orb = Int(2s+1)
    index = Dict(t => k for (k,t) in enumerate(to))
    Lmin = zeros(length(to), length(from))
    for (col, t) in enumerate(from)
        for p = 1:3
            i = t[p]
            # moving i -> i+1 never reorders the triple, so there is no fermionic sign
            if i == n_orb || (i+1) in t
                continue
            end
            m = s-i+1
            new_t = Base.setindex(t, i+1, p)
            Lmin[index[new_t], col] += sqrt(s*(s+1) - m*(m-1))
        end
    end
    return Lmin
end

function get_threebody(n_orb::Int, v_list::Vector{Int}, c_list::Vector{Float64}; tol::Float64=1e-12)
    threeBody = Vector{Tuple{ComplexF64, Tuple{Int64,Int64,Int64,Int64,Int64,Int64}}}()
    s = (n_orb-1)/2
    sectors = triple_sectors(n_orb)
    found = falses(length(v_list))
    for idx_sum in sort(collect(keys(sectors)))
        basis = sectors[idx_sum]
        M = 3s - idx_sum + 3
        # L^2 = L+L- + Lz^2 - Lz within this Lz sector
        if haskey(sectors, idx_sum+1)
            Lmin = lowering_matrix(s, basis, sectors[idx_sum+1])
            L2 = Lmin' * Lmin
        else
            L2 = zeros(length(basis), length(basis))
        end
        L2 += (M^2 - M) * I
        λ, U = eigen(Symmetric(L2))
        H = zeros(length(basis), length(basis))
        for (iv, (v, c)) in enumerate(zip(v_list, c_list))
            L = 3s - v
            sel = findall(x -> abs(x - L*(L+1)) < 1e-8 * max(1, L*(L+1)), λ)
            isempty(sel) && continue
            found[iv] = true
            V = U[:, sel]
            H += c * (V * V')
        end
        for a in eachindex(basis), b in eachindex(basis)
            if abs(H[a,b]) > tol
                push!(threeBody, (H[a,b], (basis[a]..., basis[b]...)))
            end
        end
    end
    for (iv, v) in enumerate(v_list)
        found[iv] || @warn "No three-fermion state with relative angular momentum m=$v (L=$(3s-v)) for n_orb=$n_orb"
    end
    return threeBody
end

export get_threebody
end
