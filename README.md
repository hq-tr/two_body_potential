# One-body and Two-body Potential Exact Diagonalization (ED)
by Ha Quang Trung and Fong Wan Hang

## DESCRIPTION

- This is a routine for various calculation regarding one- and two-body pseudopotentials on Landau levels. It can evaluate variational energy (with `expectvalue.jl`) and ED (with `sphere_FQHE.jl`).

- This incorporate the routines in `/sphere_FQHE_parallel/` and `/expectvalue/` into a single Julia library. It has been (mostly) optimized and is on par with these C++ routines. The Julia routines also has additional features not present in their C++ counterparts.


## USAGE
### General Instructions
- This directory contains three scripts, named "sphere_FQHE.jl", "expectvalue.jl", and "expectvalue_ED.jl" after their C++ counterpart.

	- "sphere_FQHE.jl" is a routine that perform ED using the monomials as a basis. The resulting matrix is a sparse matrix.

	- "expectvalue.jl" calculates the variational energy of a given state. The result is just a number.

	- "expectvalue_ED.jl" performs ED using a collection of states input by user as a basis. The resulting matrix is a dense matrix.

- "expectvalue.jl" is the same as "expectvalue_ED.jl" if the Hilbert space has dimension 1. However, "expectvalue.jl" contains a slightly more optimize algorithm for variational energy calculation.

- User must have the following installed:
	- The Julia programming language
	- Appropriate Julia packages
	- The QHE_Julia library

- The specific instruction for each of the routines is as follows.

### Using sphere_FQHE.jl
The general command is

`julia sphere_FQHE.jl <parameter>`

To list out all possible parameters, run help with

`julia sphere_FQHE.jl -h`

1. Specifying the basis input: There are two ways to specify a basis used to construct the Hamiltonian matrix:
	- Specifying the file name with `-f` tag. This should point to a vector file whose basis is used.
	- Specifying the number of electrons `-e` and number of orbital `-o`. Additionally, the L_z sector may be specified with the `--Lz` tag.
	Note that if both options are specified, the first option will override the second.
2. Specifying the two-body interaction: There are two ways to specify the interactions:	
	- Specifying the file name with `-i` tag. This should point to a file containing the pseudopotential description. The file should contain two columns, the first is an integer $m$ specifying the interaction $V_m$, and the second is a real number specifying the corresponding coefficient.
	- If the file name is left blank, the program will prompt the user to input the interaction data from keyboard.
	- To run the routine with no two-body interaction (one-body only), use `-i none`. Note that this means the interaction file name cannot be "none".
3. Specifying the one-body potential pin data: For now there is only one way to do this:
	- Specifying the file name with the `-p` tag. This should point to a file containing potential pin data. Each line in this file is one potential pin and should contain three columns: the first two are the azimuthal and polar angles of the pin's location, and the third is the pin's strength.
	- By default, the azimuthal and polar angles are in multiple of pi. That means specifying 0.5 in the file will mean 0.5xpi~1.571 radian. To change this, use the `--angle_multiplier` tag to change the multiplier.
4. Specifying the number of eigenvalues and eigenvectors to obtain
	- This is done by the `-n` tag. By default, this number is 5.
	
### Using expectvalue_ED.jl
(To be updated later)



## KNOWN BUGS

"expectvalue.jl" and "expectvalue_ED.jl" only works with input files with real coefficients (They do not work with complex coefficients)

