# LDLFactorizations: Factorization of Symmetric Matrices

[![Build Status](https://travis-ci.org/JuliaSmoothOptimizers/LDLFactorizations.jl.svg?branch=master)](https://travis-ci.org/JuliaSmoothOptimizers/LDLFactorizations.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/eyyrgo7qg7uxbvxm/branch/master?svg=true)](https://ci.appveyor.com/project/dpo/ldlfactorizations-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/JuliaSmoothOptimizers/LDLFactorizations.jl/badge.svg)](https://coveralls.io/github/JuliaSmoothOptimizers/LDLFactorizations.jl)

A translation of Tim Davis's Concise LDL Factorization, part of [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html).

This package is appropriate for matrices A that possess a factorization of the
form LDLᵀ without pivoting, where L is unit lower triangular and D is *diagonal* (indefinite in general), including definite and quasi-definite matrices.

LDLFactorizations.jl should not be expected to be as fast, as robust or as accurate as factorization
packages such as [HSL.jl](https://github.com/JuliaSmoothOptimizers/HSL.jl), [MUMPS.jl](https://github.com/JuliaSmoothOptimizers/MUMPS.jl) or [Pardiso.jl](https://github.com/JuliaSparse/Pardiso.jl).
Those are multifrontal and/or implement various forms of parallelism, and
employ sophisticated pivot strategies.

The main advantages of LDLFactorizations.jl are that

1. it is very short and has a small footprint;
2. it is in pure Julia, and so

   2.a. it does not require external compiled dependencies;

   2.b. it will work with multiple input data types.

Whereas MUMPS.jl, HSL.jl and Pardiso.jl only work with single and double precision
reals and complex data types, LDLFactorizations.jl accepts any numerical data type.

# Installing

```julia
julia> Pkg.add("LDLFactorizations")
```

# Usage

The only exported function is `ldl()`.
Calling `ldl()` with a dense array converts it to a sparse matrix.
A permutation ordering can be supplied: `ldl(A, p)` where `p` is an `Int`
array representing a permutation of the integers between 1 and the order
of `A`.
If no permutation is supplied, one is automatically computed using [AMD.jl](https://github.com/JuliaSmoothOptimizers/AMD.jl).
When no permutation is used, i.e., when `p = collect(1:n)`, only the upper triangle of `A` is accessed, and `ldl(A)` is the same as `ldl(triu(A))`.
Otherwise, only the upper triangle of the symmetrically-permuted `A` are accessed.
Thus, unless the identity permutation is enforced, both triangles of `A` should be supplied as input.

`ldl` returns a factorization in the form of a `LDLFactorization` object.
The `\` method is implemented for objects of type `LDLFactorization` so that
solving a linear system is as easy as
```julia
julia> LDLT = ldl(A)
julia> x = LDLT \ b  # solves Ax=b
```
The factorization can of course be reused to solve for multiple right-hand
sides.

Factors can be accessed as `LDLT.L` and `LDLT.D`, and the permutation vector as `LDLT.P`.
Because the L factor is unit lower triangular, its diagonal is not stored.
Thus the factors satisfy: PAPᵀ = (L + I) D (L + I)ᵀ.

# References

Timothy A. Davis. 2005. Algorithm 849: A concise sparse Cholesky factorization package. ACM Trans. Math. Softw. 31, 4 (December 2005), 587-591. DOI [10.1145/1114268.1114277](http://dx.doi.org/10.1145/1114268.1114277).

Like the original LDL, this package is distributed under the LGPL.

[![LGPLv3](http://www.gnu.org/graphics/lgplv3-88x31.png)](http://www.gnu.org/licenses/lgpl.html "LGPLv3")
