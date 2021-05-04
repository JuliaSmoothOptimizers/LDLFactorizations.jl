# LDLFactorizations.jl Tutorial

## A basic example

```julia
n = 10
A0 = rand(n,n)
A = A0 + A0' # A is symmetric
b = rand(n)
```

We solve the system `A x = b` using LDLFactorizations.jl :

```julia 
using LDLFactorizations, LinearAlgebra
Au = triu(A) # get upper triangle
LDL = ldl(Symmetric(Au, :U))
x = LDL \ b
```

## A more performance-focused example

We build a problem with sparse arrays.

```julia
using SparseArrays
n = 100
A0 = sprand(Float64, n, n, 0.1)
A = A0 * A0' + I # A is symmetric
b = rand(n)
```

Now if we want to use the factorization to solve multiple systems that have 
the same sparsity pattern as A, we only have to use `ldl_analyze` once.

```julia 
Au = triu(A) # get upper triangle
x = zeros(n)

LDL = ldl_analyze(Symmetric(Au, :U)) # symbolic analysis
ldl_factorize!(Symmetric(Au, :U), LDL) # factorization
ldiv!(x, LDL, b) # solve in-place (we could use ldiv!(LDL, b) if we want to overwrite b)

Au.nzval .+= 1.0 # modify Au without changing the sparsity pattern
ldl_factorize!(Symmetric(Au, :U), LDL) 
ldiv!(x, LDL, b)
```

## Dynamic Regularization

When the matrix to factorize is nearly singular and the factorization encounters zero pivots, 
if we know the signs of the pivots and if they are clustered by signs (for example, the 
`n_d` first pivots are positive and the other pivots are negative), we can use:

```julia
ϵ = sqrt(eps())
LDL = ldl_analyze(Symmetric(triu(A), :U))
LDL.tol = ϵt
LDL.n_d = 10
LDL.r1 = ϵ # if any of the n_d first pivots D[i] < ϵ, then D[i] += LDL.r1    
LDL.r2 = -ϵ # if any of the n - n_d last pivots D[i] < ϵ, then D[i] += LDL.r2 
ldl_factorize!(Symmetric(triu(A), :U), LDL)
```