using LDLFactorizations, Test, LinearAlgebra, SparseArrays

@testset "factorizable" begin
  # this matrix possesses an LDLᵀ factorization without pivoting
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

  LDLT = ldl(A)
  x = LDLT \ b

  r = A * x - b
  @test norm(r) ≤ ϵ * norm(b)

  y = collect(0.1:0.1:1)
  @test norm(x - y) ≤ ϵ * norm(y)

  x2 = copy(b)
  ldiv!(LDLT, x2)

  r2 = A * x2 - b
  @test norm(r2) ≤ ϵ * norm(b)

  @test norm(x2 - y) ≤ ϵ * norm(y)

  # test properties
  @test eltype(LDLT) == eltype(A)
  @test nnz(LDLT) == nnz(LDLT.L) + length(LDLT.d)
  @test size(LDLT) == size(A)
  @test typeof(LDLT.D) <: Diagonal
  @test propertynames(LDLT) == (:L, :D, :P)
end

@testset "not_factorizable" begin
  # this matrix does not possess an LDLᵀ factorization without pivoting
  A = [ 0 1
        1 1 ]
  @test_throws LDLFactorizations.SQDException ldl(A, [1, 2])

  A = Symmetric(sparse(triu(A)))
  S = ldl_analyze(A)
  @test_throws LDLFactorizations.SQDException ldl_factorize!(A, S)
end

@testset "sparse" begin
  for Ti in (Int32, Int), Tf in (Float32, Float64, BigFloat)
    A = sparse(Ti[1, 2, 1, 2], Ti[1, 1, 2, 2], Tf[10, 2, 2, 5])
    z = ones(Tf, 2)
    b = A * z
    LDLT = ldl(A)
    x = LDLT \ b
    @test norm(x - z) ≤ sqrt(eps(Tf)) * norm(z)
    r = A * x - b
    @test norm(r) ≤ sqrt(eps(Tf)) * norm(b)

    x2 = copy(b)
    ldiv!(LDLT, x2)
    @test norm(x2 - z) ≤ sqrt(eps(Tf)) * norm(z)
    r2 = A * x2 - b
    @test norm(r2) ≤ sqrt(eps(Tf)) * norm(b)

    y = similar(b)
    ldiv!(y, LDLT, b)
    @test norm(y - z) ≤ sqrt(eps(Tf)) * norm(z)
    r2 = A * y - b
    @test norm(r2) ≤ sqrt(eps(Tf)) * norm(b)

    @test nnz(LDLT) == nnz(LDLT.L) + length(LDLT.d)

    # test with separate analyze/factorize phases
    # A = Symmetric(sparse(triu(A)))
    S = ldl_analyze(A)
    ldl_factorize!(A, S)
    x2 = copy(b)
    ldiv!(LDLT, x2)
    @test norm(x2 - z) ≤ sqrt(eps(Tf)) * norm(z)
    r2 = A * x2 - b
    @test norm(r2) ≤ sqrt(eps(Tf)) * norm(b)
  end
end

@testset "factorizable_upper" begin
  # Using only the upper triangle tests
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
  LDLT = ldl(A)

  A_upper = Symmetric(triu(A), :U)
  LDLT_upper = ldl(A_upper)
  x = LDLT_upper \ b

  y = collect(0.1:0.1:1)
  @test norm(x - y) ≤ ϵ * norm(y)

  r = A * x - b
  @test norm(r) ≤ ϵ * norm(b)

  x2 = copy(b)
  ldiv!(LDLT, x2)
  @test norm(x2 - y) ≤ ϵ * norm(y)

  r2 = A * x2 - b
  @test norm(r2) ≤ ϵ * norm(b)

  @test nnz(LDLT_upper) == nnz(LDLT_upper.L) + length(LDLT_upper.d)

  # test with separate analyze/factorize phases
  # A = Symmetric(sparse(triu(A)))
  S = ldl_analyze(A)
  ldl_factorize!(A, S)
  x2 = copy(b)
  ldiv!(LDLT, x2)
  @test norm(x2 - y) ≤ ϵ * norm(y)
  r2 = A * x2 - b
  @test norm(r2) ≤ ϵ * norm(b)

  # test with a permutation of a different int type
  p = Int32.(collect(size(A,1):-1:1))
  LDLT_upper = ldl(A_upper, p)
  x = LDLT_upper \ b

  y = collect(0.1:0.1:1)
  @test norm(x - y) ≤ ϵ * norm(y)

  r = A * x - b
  @test norm(r) ≤ ϵ * norm(b)

  x2 = copy(b)
  ldiv!(LDLT, x2)
  @test norm(x2 - y) ≤ ϵ * norm(y)

  r2 = A * x2 - b
  @test norm(r2) ≤ ϵ * norm(b)

  # Tests with multiple right-hand sides
  B = 1.0 * [ i + j for j = 1 : 10, i = 0 : 3]
  X = A \ B
  Y = similar(B)
  ldiv!(Y, LDLT, B)
  @test norm(Y - X) ≤ ϵ * norm(X)
end

@testset "not_factorizable_upper" begin
  # this matrix does not possess an LDLᵀ factorization without pivoting
  A = triu([ 0 1
            1 1 ])
  @test_throws LDLFactorizations.SQDException ldl(A, [1, 2])
end

@testset "sparse_upper" begin
  for Ti in (Int32, Int), Tf in (Float32, Float64, BigFloat)
    A = sparse(Ti[1, 2, 1, 2], Ti[1, 1, 2, 2], Tf[10, 2, 2, 5])
    A_upper = Symmetric(triu(A), :U)
    b = A * ones(Tf, 2)
    LDLT = ldl(A_upper)
    x = LDLT \ b
    r = A * x - b
    @test norm(r) ≤ sqrt(eps(Tf)) * norm(b)

    x2 = copy(b)
    ldiv!(LDLT, x2)
    r2 = A * x2 - b
    @test norm(r2) ≤ sqrt(eps(Tf)) * norm(b)

    @test nnz(LDLT) == nnz(LDLT.L) + length(LDLT.d)
  end
end

@testset "positive_semidefinite" begin
  A = [0.   0.   0.   0.   0.   0.   0.   0.   4.   0.
       0.   0.   0.   0.   0.   0.   0.   0.   5.   0.
       2.   4.   5.   -2   4.   1.   2.   2.   2.   0.
       0.   0.   0.   0.   1.   9.   9.   1.   7.   1.
       0.   0.   0.   0.   0.   0.   0.   0.   1.   0.
       1.   3.   2.   1.   4.   3.   1.   0.   0.   7.
       -3.  8.   0.   0.   0.   0.   -2.  0.   0.   1.
       0.   0.   0.   5.   7.   9.   0.   2.   7.   1.
       3.   2.   0.   0.   0.   0.   1.   3.   3.   2.
       0.   0.   0.   0.  -3   -4    0.   0.   0.   0. ]
  M = A * A'  # det(A) = 0 => M positive semidefinite
  b = M * ones(10)
  x = copy(b)
  S = ldl_analyze(Symmetric(triu(M), :U))
  S = ldl_factorize!(Symmetric(triu(M), :U), S, tol=1e-8, r1=0., r2=1e-8, n_d=0)
  x = ldiv!(S, x)
  r = M * x - b
  @test norm(r) ≤ sqrt(eps()) * norm(b)
end

@testset "SQD" begin
  A = [0.   0.   0.   0.   0.   0.   0.   0.   4.   0.
       0.   0.   0.   0.   0.   0.   0.   0.   5.   0.
       2.   4.   5.   -2   4.   1.   2.   2.   2.   0.
       0.   0.   0.   0.   1.   9.   9.   1.   7.   1.
       0.   0.   0.   0.   0.   0.   0.   0.   1.   0.
       1.   3.   2.   1.   4.   3.   1.   0.   0.   7.
       -3.  8.   0.   0.   0.   0.   -2.  0.   0.   1.
       0.   0.   0.   5.   7.   9.   0.   2.   7.   1.
       3.   2.   0.   0.   0.   0.   1.   3.   3.   2.
       0.   0.   0.   0.  -3   -4    0.   0.   0.   0. ]
  M = spzeros(20, 20)
  M[1:10, 1:10] = -A * A'
  M[11:20, 11:20] = A * A'
  # M = [-A*A'    0
  #        0     A*A'] where A*A' is symmetric positive semidefinite
  b = M * ones(20)
  x = copy(b)
  S = ldl_analyze(Symmetric(triu(M), :U))
  S = ldl_factorize!(Symmetric(triu(M), :U), S, tol=1e-8, r1=-1e-8, r2=1e-8, n_d=0)
  x = ldiv!(S, x)
  r = M * x - b
  @test norm(r) ≤ sqrt(eps()) * norm(b)
end
