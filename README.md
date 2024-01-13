This is a routine for various calculation regarding two-body pseudopotentials. It can evaluate variational energy (with `expectvalue.jl`) and ED (with `sphere_FQHE.jl`).
This incorporate the routines in `/sphere_FQHE_parallel/` and `/expectvalue/` into a single Julia library. However, at the moment the Julia routine is far from optimized and can't stand up to the other two C++ routines.

The ultimate goal is to incorporate two-body and one-body potentials into a single ED routine.
