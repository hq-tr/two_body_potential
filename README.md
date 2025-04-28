=============================================================================


            	Two-body potential exact diagonalization (ED)
            		by Ha Quang Trung and Fong Wan Hang


=============================================================================

===== DESCRIPTION
• This is a routine for various calculation regarding two-body pseudopotentials. It can evaluate variational energy (with `expectvalue.jl`) and ED (with `sphere_FQHE.jl`).

• This incorporate the routines in `/sphere_FQHE_parallel/` and `/expectvalue/` into a single Julia library. However, at the moment the Julia routine is far from optimized and can't stand up to the other two C++ routines.

• The ultimate goal is to incorporate two-body and one-body potentials into a single ED routine.


===== USAGE

• This directory contains three scripts, named "sphere_FQHE.jl", "expectvalue.jl", and "expectvalue_ED.jl" after their C++ counterpart.

	• "sphere_FQHE.jl" is a routine that perform ED using the monomials as a basis. The resulting matrix is a sparse matrix.

	• "expectvalue.jl" calculates the variational energy of a given state. The result is just a number.

	• "expectvalue_ED.jl" performs ED using a collection of states input by user as a basis. The resulting matrix is a dense matrix.

• "expectvalue.jl" is the same as "expectvalue_ED.jl" if the Hilbert space has dimension 1. However, "expectvalue.jl" contains a slightly more optimize algorithm for variational energy calculation.


===== KNOWN BUGS

"expectvalue.jl" and "expectvalue_ED.jl" only works with input files with real coefficients (They do not work with complex coefficients)

