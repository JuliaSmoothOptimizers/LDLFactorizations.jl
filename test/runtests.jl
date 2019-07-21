using LDLFactorizations, Test, LinearAlgebra, SparseArrays

# this matrix possesses an LDL factorization without pivoting
A = [ 1.7     0     0     0     0     0     0     0   .13     0
        0    1.     0     0   .02     0     0     0     0   .01
        0     0   1.5     0     0     0     0     0     0     0
        0     0     0   1.1     0     0     0     0     0     0
        0   .02     0     0   2.6     0   .16   .09   .52   .53
        0     0     0     0     0   1.2     0     0     0     0
        0     0     0     0   .16     0   1.3     0     0   .56
        0     0     0     0   .09     0     0   1.6   .11     0
      .13     0     0     0   .52     0     0   .11   1.4     0
        0   .01     0     0   .53     0   .56     0     0   3.1 ]
b = [.287, .22, .45, .44, 2.486, .72, 1.55, 1.424, 1.621, 3.759]
ϵ = sqrt(eps(eltype(A)))

y = collect(0.1:0.1:1)

LDLT = ldl(A)
x = LDLT \ b
r = A * x - b
@test norm(r) ≤ ϵ * norm(b)
@test norm(x - y) ≤ ϵ * norm(y)

# test with input containing only upper triangular part
LDLT = ldl(sparse(UpperTriangular(A)), onlyupper = true)
x = LDLT \ b
r = A * x - b
@test norm(r) ≤ ϵ * norm(b)
@test norm(x - y) ≤ ϵ * norm(y)

# this matrix does not possess an LDL factorization without pivoting
A = [ 0 1
      1 1 ]
@test_throws LDLFactorizations.SQDException ldl(A, [1, 2])

# Sparse tests
for Ti in (Int32, Int), Tf in (Float32, Float64, BigFloat)
  A = sparse(Ti[1, 2, 1, 2], Ti[1, 1, 2, 2], Tf[10, 2, 2, 5])
  b = A * ones(Tf, 2)
  P = Ti[2; 1]

  for LDLT in [ldl(A), ldl(A, P), ldl(sparse(UpperTriangular(A)), P, onlyupper = true)]
    x = LDLT \ b
    r = A * x - b
    @test norm(r) ≤ sqrt(eps(Tf)) * norm(b)
  end
end
