module LDLFactorizations

export ldl, ldl_analyze, ldl_factorize!, \, ldiv!, nnz

using AMD, LinearAlgebra, SparseArrays

mutable struct SQDException <: Exception
  msg::String
end

"""
    col_symb!(n, Ap, Ai, Cp, w, Pinv)
Compute the sparse structure of missing elements of the upper triangle of PAPt. Nonzero elements have to verify Pinv[i] < Pinv[j] where i is the
row index and j the column index. Those elements are the nonzeros of the lower triangle of A that will be in the upper triangle of PAPt (after permutation)
# Arguments
- `n::Ti`: number of columns of the matrix
- `Ap::Vector{Ti}`: colptr of the matrix to factorize (CSC format)
- `Ai::Vector{Ti}`: rowval of the matrix to factorize (CSC format)
- `Cp::Vector{Ti}`: colptr of the lower triangle (to be modified)
- `w::Vector{Ti}`: work array
- `Pinv::Vector{Ti}`: inverse permutation of P. PAPt is the matrix to factorize (CSC format)
"""
function col_symb!(n, Ap, Ai, Cp, w, Pinv)
  fill!(w, 0)
  @inbounds for j = 1:n
    @inbounds for p = Ap[j] : (Ap[j+1]-1)
      i = Ai[p]
      i >= j && break  # only upper part
      Pinv[i] < Pinv[j] && continue  # store only what will be used during the factorization
      w[i] += 1  # count entries
    end
  end

  Cp[1] = 1
  @inbounds for i = 1:n  # cumulative sum
    Cp[i+1] = w[i] + Cp[i]
    w[i] = Cp[i]
  end
end

"""
    col_num!(n, Ap, Ai, Ci, w, Pinv)
    
Compute the rowval and values of missing elements of the upper triangle of PAPt. Nonzero elements have to verify Pinv[i] ≥ Pinv[j] where i is the
row index and j the column index. Those elements are the nonzeros of the lower triangle of A that will be in the upper triangle of PAPt (after permutation)
# Arguments
- `n::Ti`: number of columns of the matrix
- `Ap::Vector{Ti}`: colptr of the matrix to factorize (CSC format)
- `Ai::Vector{Ti}`: rowval of the matrix to factorize (CSC format)
- `Ci::Vector{Ti}`: rowval of the lower triangle
- `w::Vector{Ti}`: work array
- `Pinv::Vector{Ti}`: inverse permutation of P. PAPt is the matrix to factorize (CSC format)
"""
function col_num!(n, Ap, Ai, Ci, w, Pinv)
  @inbounds for j = 1:n
    @inbounds for p = Ap[j] : (Ap[j+1]-1)
      i = Ai[p]
      i >= j  && break  # only upper part
      Pinv[i] < Pinv[j] && continue  # store only what will be used during the factorization
      Ci[w[i]] = j
      w[i] += 1
    end
  end
end

function ldl_symbolic_upper!(n, Ap, Ai, Cp, Ci , Lp, parent, Lnz, flag, P, Pinv)
  @inbounds for k = 1:n
    parent[k] = -1
    flag[k] = k
    Lnz[k] = 0
    pk = P[k]
    @inbounds for p = Ap[pk] : (Ap[pk+1] - 1)
      i = Pinv[Ai[p]]
      i ≥ k && continue
      @inbounds while flag[i] != k
        if parent[i] == -1
          parent[i] = k
        end
        Lnz[i] += 1
        flag[i] = k
        i = parent[i]
      end
    end

    # Missing nonzero elements of the upper triangle
    @inbounds for ind = Cp[pk] : (Cp[pk+1] - 1)
      i = Pinv[Ci[ind]]
      i > k && continue
      @inbounds while flag[i] != k
        if parent[i] == -1
          parent[i] = k
        end
      Lnz[i] += 1
      flag[i] = k
      i = parent[i]
      end
    end
  end
  Lp[1] = 1
  @inbounds for k = 1:n
    Lp[k+1] = Lp[k] + Lnz[k]
  end
end

function ldl_symbolic!(n, Ap, Ai, Lp, parent, Lnz, flag, P, Pinv)
  @inbounds for k = 1:n
  parent[k] = -1
  flag[k] = k
  Lnz[k] = 0
  pk = P[k]
  @inbounds for p = Ap[pk] : (Ap[pk+1] - 1)
    i = Pinv[Ai[p]]
    i ≥ k && continue
      @inbounds while flag[i] != k
        if parent[i] == -1
          parent[i] = k
        end
        Lnz[i] += 1
        flag[i] = k
        i = parent[i]
        end
      end
    end
    Lp[1] = 1
    @inbounds for k = 1:n
    Lp[k+1] = Lp[k] + Lnz[k]
  end
end

function ldl_numeric_upper!(n, Ap, Ai, Ax, Cp, Ci, Lp, parent, Lnz, Li, Lx, D, Y,
                            pattern, flag, P, Pinv; dynamic_args...)
  T = eltype(Ax)
  if length(dynamic_args) != 0
    dynamic_regul = true
    dynamic_dict = Dict(dynamic_args)
  else
    dynamic_regul = false
  end
  @inbounds for k = 1:n
    Y[k] = 0
    top = n+1
    flag[k] = k
    Lnz[k] = 0
    pk = P[k]
    @inbounds for p = Ap[pk] : (Ap[pk+1] - 1)
      i = Pinv[Ai[p]]
      i > k && continue
      Y[i] += Ax[p]
      len = 1
      @inbounds while flag[i] != k
        pattern[len] = i
        len += 1
        flag[i] = k
        i = parent[i]
      end
      @inbounds while len > 1
        top -= 1
        len -= 1
        pattern[top] = pattern[len]
      end
    end
    # missing non zero elements of the upper triangle
    @inbounds for ind = Cp[pk] : (Cp[pk+1] - 1)
      i2 = Ci[ind]
      i = Pinv[i2]
      i > k && continue
      @inbounds for p = Ap[i2] : (Ap[i2+1] - 1)
        Ai[p] < pk && continue
        Y[i] += Ax[p]
        len = 1
        @inbounds while flag[i] != k
          pattern[len] = i
          len += 1
          flag[i] = k
          i = parent[i]
        end
        @inbounds while len > 1
          top -= 1
          len -= 1
          pattern[top] = pattern[len]
        end
        break
      end
    end
    D[k] = Y[k]
    Y[k] = 0
    @inbounds while top ≤ n
      i = pattern[top]
      yi = Y[i]
      Y[i] = 0
      @inbounds for p = Lp[i] : (Lp[i] + Lnz[i] - 1)
        Y[Li[p]] -= Lx[p] * yi
      end
      p = Lp[i] + Lnz[i]
      l_ki = yi / D[i]
      D[k] -= l_ki * yi
      Li[p] = k
      Lx[p] = l_ki
      Lnz[i] += 1
      top += 1
    end
    if dynamic_regul && abs(D[k]) < eps(T) * dynamic_dict[:Amax]
      if P[k] <= dynamic_dict[:n_d]
        D[k] += dynamic_dict[:r1]
      else
        D[k] += dynamic_dict[:r2]
      end
    end
    D[k] == 0 && throw(SQDException("matrix does not possess a LDL' factorization for this permutation"))
  end
end

function ldl_numeric!(n, Ap, Ai, Ax, Lp, parent, Lnz, Li, Lx, D, Y,
                      pattern, flag, P, Pinv)
  @inbounds for k = 1:n
    Y[k] = 0
    top = n + 1
    flag[k] = k
    Lnz[k] = 0
    pk = P[k]
    @inbounds for p = Ap[pk] : (Ap[pk+1] - 1)
      i = Pinv[Ai[p]]
      i > k && continue
      Y[i] += Ax[p]
      len = 1
      @inbounds  while flag[i] != k
        pattern[len] = i
        len += 1
        flag[i] = k
        i = parent[i]
      end
      @inbounds while len > 1
        top -= 1
        len -= 1
        pattern[top] = pattern[len]
      end
    end
    D[k] = Y[k]
    Y[k] = 0
    @inbounds while top ≤ n
      i = pattern[top]
      yi = Y[i]
      Y[i] = 0
      @inbounds for p = Lp[i] : (Lp[i] + Lnz[i] - 1)
        Y[Li[p]] -= Lx[p] * yi
      end
      p = Lp[i] + Lnz[i]
      l_ki = yi / D[i]
      D[k] -= l_ki * yi
      Li[p] = k
      Lx[p] = l_ki
      Lnz[i] += 1
      top += 1
    end
    D[k] == 0 && throw(SQDException("matrix does not possess a LDL' factorization for this permutation"))
  end
end

# solve functions for a single rhs
function ldl_lsolve!(n, x::AbstractVector, Lp, Li, Lx)
  @inbounds for j = 1:n
    xj = x[j]
    @inbounds for p = Lp[j] : (Lp[j+1] - 1)
      x[Li[p]] -= Lx[p] * xj
    end
  end
  return x
end

function ldl_dsolve!(n, x::AbstractVector, D)
  @inbounds for j = 1:n
    x[j] /= D[j]
  end
  return x
end

function ldl_ltsolve!(n, x::AbstractVector, Lp, Li, Lx)
  @inbounds for j = n:-1:1
    xj = x[j]
    @inbounds for p = Lp[j] : (Lp[j+1] - 1)
      xj -= Lx[p] * x[Li[p]]
    end
    x[j] = xj
  end
  return x
end

function ldl_solve!(n, b::AbstractVector, Lp, Li, Lx, D, P)
  @views y = b[P]
  ldl_lsolve!(n, y, Lp, Li, Lx)
  ldl_dsolve!(n, y, D)
  ldl_ltsolve!(n, y, Lp, Li, Lx)
  return b
end

# solve functions for multiple rhs
function ldl_lsolve!(n, X::AbstractMatrix{T}, Lp, Li, Lx) where T
  @inbounds for j = 1:n
    @inbounds for p = Lp[j] : (Lp[j+1] - 1)
      for k ∈ axes(X, 2)
        X[Li[p], k] -= Lx[p] * X[j, k]
      end
    end
  end
  return X
end

function ldl_dsolve!(n, X::AbstractMatrix{T}, D) where T
  @inbounds for j = 1:n
    for k ∈ axes(X, 2)
      X[j, k] /= D[j]
    end
  end
  return X
end

function ldl_ltsolve!(n, X::AbstractMatrix{T}, Lp, Li, Lx) where T
  @inbounds for j = n:-1:1
    @inbounds for p = Lp[j] : (Lp[j+1] - 1)
      for k ∈ axes(X, 2)
        X[j, k] -= Lx[p] * X[Li[p], k]
      end
    end
  end
  return X
end

function ldl_solve!(n, B::AbstractMatrix{T}, Lp, Li, Lx, D, P) where T
  @views Y = B[P, :]
  ldl_lsolve!(n, Y, Lp, Li, Lx)
  ldl_dsolve!(n, Y, D)
  ldl_ltsolve!(n, Y, Lp, Li, Lx)
  return B
end

# a simplistic type for LDLᵀ factorizations so we can do \ and separate analyze/factorize
mutable struct LDLFactorization{T<:Real,Ti<:Integer,Tn<:Integer,Tp<:Integer}
  __analyzed::Bool
  __factorized::Bool
  __upper::Bool
  n::Tn
  # fields related to symbolic analysis
  parent::Vector{Ti}
  Lnz::Vector{Ti}
  flag::Vector{Ti}
  P::Vector{Tp}
  pinv::Vector{Tp}
  Lp::Vector{Ti}
  Cp::Vector{Ti}
  Ci::Vector{Ti}
  # fields related to numerical factorization
  Li::Vector{Ti}
  Lx::Vector{T}
  d::Vector{T}
  Y::Vector{T}
  pattern::Vector{Ti}
end

# perform symbolic analysis so it can be reused
function ldl_analyze(A::Symmetric{T,SparseMatrixCSC{T,Ti}}, P::Vector{Tp}) where {T<:Real,Ti<:Integer,Tp<:Integer}
  A.uplo == 'U' || error("upper triangle must be supplied")
  n = size(A, 1)
  n == size(A, 2) || throw(DimensionMismatch("matrix must be square"))
  n == length(P) || throw(DimensionMismatch("permutation size mismatch"))

  # allocate space for symbolic analysis
  parent = Vector{Ti}(undef, n)
  Lnz = Vector{Ti}(undef, n)
  flag = Vector{Ti}(undef, n)
  pinv = Vector{Tp}(undef, n)
  Lp = Vector{Ti}(undef, n+1)

  # Compute inverse permutation
  @inbounds for k = 1:n
    pinv[P[k]] = k
  end

  Cp = Vector{Ti}(undef, n + 1)
  col_symb!(n, A.data.colptr, A.data.rowval, Cp, Lp, pinv)
  Ci = Vector{Ti}(undef, Cp[end] - 1)
  col_num!(n, A.data.colptr, A.data.rowval, Ci, Lp, pinv)

  # perform symbolic analysis
  ldl_symbolic_upper!(n, A.data.colptr, A.data.rowval, Cp, Ci, Lp, parent, Lnz, flag, P, pinv)

  # space for numerical factorization will be allocated later
  Li = Ti[]
  Lx = T[]
  D = T[]
  Y = T[]
  pattern = Ti[]
  return LDLFactorization(true, false, true, n, parent, Lnz, flag, P, pinv, Lp, Cp, Ci, Li, Lx, D, Y, pattern)
end

# convert dense to sparse
ldl_analyze(A::Symmetric{T,Array{T,2}}) where T<:Real = ldl_analyze(Symmetric(sparse(A.data)))
ldl_analyze(A::Symmetric{T,Array{T,2}}, P) where T<:Real = ldl_analyze(Symmetric(sparse(A.data)), P)

# use AMD permuation by default
ldl_analyze(A::Symmetric{T,SparseMatrixCSC{T,Ti}}) where {T<:Real,Ti<:Integer} = ldl_analyze(A, amd(A))

function ldl_factorize!(A::Symmetric{T,SparseMatrixCSC{T,Ti}},
                        S::LDLFactorization{T,Ti,Tn,Tp};
                        dynamic_args...) where {T<:Real,Ti<:Integer,Tn<:Integer,Tp<:Integer}
  # dynamic_args... is used to perform dynamic regularization of S.d, for example if T=Float64:
  # ldl_factorize!(A, S, Amax=100.0, n_d=20, r1=-eps()^(1/4), r2=eps()^(1/2))
  # in this case, if i<=n_d and abs(S.d[i]) <= Amax*eps() (without permuting), S.d[i]+=r1
  # (resp. S.d[i]+=r2 if i>n_d)
  S.__analyzed || error("perform symbolic analysis prior to numerical factorization")
  n = size(A, 1)
  n == S.n || throw(DimensionMismatch("matrix size is inconsistent with symbolic analysis object"))

  # allocate space for numerical factorization if not already done
  if !(S.__factorized)
    S.Li = Vector{Ti}(undef, S.Lp[n] - 1)
    S.Lx = Vector{T}(undef, S.Lp[n] - 1)
    S.d = Vector{T}(undef, n)
    S.Y = Vector{T}(undef, n)
    S.pattern = Vector{Ti}(undef, n)
  end

  # perform numerical factorization
  ldl_numeric_upper!(S.n, A.data.colptr, A.data.rowval, A.data.nzval,
                     S.Cp, S.Ci, S.Lp, S.parent, S.Lnz, S.Li, S.Lx, S.d, S.Y, S.pattern, S.flag, S.P, S.pinv;
                     dynamic_args...)
  return S
end

# convert dense to sparse
ldl_factorize!(A::Symmetric{T,Array{T,2}}, S::LDLFactorization) where T<:Real = ldl_factorize!(Symmetric(sparse(A.data)), S)

# symmetric matrix input
function ldl(sA::Symmetric{T,SparseMatrixCSC{T,Ti}}, P::Vector{Tp}) where {T<:Real,Ti<:Integer,Tp<:Integer}
  sA.uplo == 'L' && error("matrix must contain the upper triangle")
  # ldl(sA.data, args...; upper = true )
  S = ldl_analyze(sA, P)
  ldl_factorize!(sA, S)
end

ldl(sA::Symmetric{T,Array{T,2}}) where T<:Real = ldl(Symmetric(sparse(sA.data)))
ldl(sA::Symmetric{T,Array{T,2}}, P) where T<:Real = ldl(Symmetric(sparse(sA.data)), P)
ldl(sA::Symmetric{T,SparseMatrixCSC{T,Ti}}) where {T<:Real,Ti<:Integer} = ldl(sA, amd(sA))

# convert dense to sparse
ldl(A::Array{T,2}) where T<:Real = ldl(sparse(A))
ldl(A::Array{T,2}, P) where T<:Real = ldl(sparse(A), P)

# use AMD permutation by default
ldl(A::SparseMatrixCSC{T,Ti}) where {T<:Real,Ti<:Integer} = ldl(A, amd(A))

# use ldl(A, collect(1:n)) to suppress permutation
function ldl_analyze(A::SparseMatrixCSC{T,Ti}, P::Vector{Tp}) where {T<:Real,Ti<:Integer,Tp<:Integer}
  n = size(A, 1)
  n == size(A, 2) || throw(DimensionMismatch("matrix must be square"))
  n == length(P) || throw(DimensionMismatch("permutation size mismatch"))

  # allocate space for symbolic analysis
  parent = Vector{Ti}(undef, n)
  Lnz = Vector{Ti}(undef, n)
  flag = Vector{Ti}(undef, n)
  pinv = Vector{Tp}(undef, n)
  Lp = Vector{Ti}(undef, n+1)
  Cp = Ti[]
  Ci = Ti[]

  # Compute inverse permutation
  @inbounds for k = 1:n
    pinv[P[k]] = k
  end

  ldl_symbolic!(n, A.colptr, A.rowval, Lp, parent, Lnz, flag, P, pinv)

  # space for numerical factorization will be allocated later
  Li = Ti[]
  Lx = T[]
  D = T[]
  Y = T[]
  pattern = Ti[]
  return LDLFactorization(true, false, false, n, parent, Lnz, flag, P, pinv, Lp, Cp, Ci, Li, Lx, D, Y, pattern)
end

# convert dense to sparse
ldl_analyze(A::Array{T,2}) where T<:Real = ldl_analyze(sparse(A))
ldl_analyze(A::Array{T,2}, P) where T<:Real = ldl_analyze(sparse(A), P)

# use AMD permuation by default
ldl_analyze(A::SparseMatrixCSC{T,Ti}) where {T<:Real,Ti<:Integer} = ldl_analyze(A, amd(A))

function ldl_factorize!(A::SparseMatrixCSC{T,Ti},
                        S::LDLFactorization{T,Ti,Tn,Tp}) where {T<:Real,Ti<:Integer,Tn<:Integer,Tp<:Integer}
  S.__upper && error("symbolic analysis was performed for a Symmetric{} matrix")
  S.__analyzed || error("perform symbolic analysis prior to numerical factorization")
  n = size(A, 1)
  n == S.n || throw(DimensionMismatch("matrix size is inconsistent with symbolic analysis object"))

  if !(S.__factorized)
    S.Li = Vector{Ti}(undef, S.Lp[n] - 1)
    S.Lx = Vector{T}(undef, S.Lp[n] - 1)
    S.d = Vector{T}(undef, n)
    S.Y = Vector{T}(undef, n)
    S.pattern = Vector{Ti}(undef, n)
  end

  ldl_numeric!(S.n, A.colptr, A.rowval, A.nzval, S.Lp, S.parent, S.Lnz,
               S.Li, S.Lx, S.d, S.Y, S.pattern, S.flag, S.P, S.pinv)
  return S
end

# convert dense to sparse
ldl_factorize!(A::Array{T,2}, S::LDLFactorization) where T<:Real = ldl_factorize!(sparse(A), S)

function ldl(A::SparseMatrixCSC, P::Vector{Tp}) where Tp <: Integer
  S = ldl_analyze(A, P)
  ldl_factorize!(A, S)
end

import Base.(\)
function (\)(LDL::LDLFactorization{T,Ti,Tn,Tp}, b::AbstractVector{T}) where {T<:Real,Ti<:Integer,Tn<:Integer,Tp<:Integer}
  y = copy(b)
  ldl_solve!(LDL.n, y, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)
end

function (\)(LDL::LDLFactorization{T,Ti,Tn,Tp}, B::AbstractMatrix{T}) where {T<:Real,Ti<:Integer,Tn<:Integer,Tp<:Integer}
  Y = copy(B)
  ldl_solve!(LDL.n, Y, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)
end

import LinearAlgebra.ldiv!
@inline ldiv!(LDL::LDLFactorization{T,Ti,Tn,Tp}, b::AbstractVector{T}) where {T<:Real,Ti<:Integer,Tn<:Integer,Tp<:Integer} =
  ldl_solve!(LDL.n, b, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)

@inline ldiv!(LDL::LDLFactorization{T,Ti,Tn,Tp}, B::AbstractMatrix{T}) where {T<:Real,Ti<:Integer,Tn<:Integer,Tp<:Integer} =
  ldl_solve!(LDL.n, B, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)

function ldiv!(y::AbstractVector{T}, LDL::LDLFactorization{T,Ti,Tn,Tp}, b::AbstractVector{T}) where {T<:Real,Ti<:Integer,Tn<:Integer,Tp<:Integer}
  y .= b
  ldl_solve!(LDL.n, y, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)
end

function ldiv!(Y::AbstractMatrix{T}, LDL::LDLFactorization{T,Ti}, B::AbstractMatrix{T}) where {T<:Real,Ti<:Integer,Tn<:Integer,Tp<:Integer}
  Y .= B
  ldl_solve!(LDL.n, Y, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)
end

Base.eltype(LDL::LDLFactorization) = eltype(LDL.d)
Base.size(LDL::LDLFactorization) = (LDL.n, LDL.n)
SparseArrays.nnz(LDL::LDLFactorization) = length(LDL.Lx) + length(LDL.d)

# user-friendly factors
@inline function Base.getproperty(LDL::LDLFactorization, prop::Symbol)
  if prop == :L
    # TODO: permute and return UnitLowerTriangular()
    return SparseMatrixCSC(LDL.n, LDL.n, LDL.Lp, LDL.Li, LDL.Lx)
  elseif prop == :D
    return Diagonal(LDL.d)
  else
    getfield(LDL, prop)
  end
end

Base.propertynames(LDL::LDLFactorization) = (:L, :D, :P)

end  # module
