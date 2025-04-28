include("/home/trung/_qhe-julia/HilbertSpace.jl")
include("/home/trung/_qhe-julia/FQH_state_v2.jl")
using .HilbertSpaceGenerator
using .FQH_states

println("How many electrons? ")
Ne = parse(Int,readline())
println("How many orbitals? ")
No = parse(Int,readline())
println("Input L_z sector on the sphere (leave blank for all L_z sectors.")
Lz_text = readline()

if length(Lz_text) > 0
	L_z = parse(Int, Lz_text)
	basis = fullhilbertspace(Ne,No,L_z)
else
	basis = fullhilbertspace(Ne,No)
end

dim = length(basis)
println("The dimension is $(dim).")

coef = zeros(dim)

if length(Lz_text) > 0
	printwf(FQH_state(basis,coef);fname="b_$(Ne)e_$(No)o_Lz_$(L_z)")
else
	printwf(FQH_state(basis,coef);fname="b_$(Ne)e_$(No)o")
end

