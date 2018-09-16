module LDLFactorizations

export ldl, \

using AMD, SparseArrays

mutable struct SQDException <: Exception
  msg::String
end

function ldl_symbolic!(n, Ap, Ai, Lp, parent, Lnz, flag, P, Pinv)
  # compute inverse permutation
  @inbounds for k = 1:n
    Pinv[P[k]] = k
  end

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

function ldl_numeric!(n, Ap, Ai, Ax, Lp, parent, Lnz, Li, Lx, D, Y,
                      pattern, flag, P, Pinv)
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

function ldl_lsolve!(n, x, Lp, Li, Lx)
  @inbounds for j = 1:n
    xj = x[j]
    @inbounds for p = Lp[j] : (Lp[j+1] - 1)
      x[Li[p]] -= Lx[p] * xj
    end
  end
  return x
end

function ldl_dsolve!(n, x, D)
  @inbounds for j = 1:n
    x[j] /= D[j]
  end
  return x
end

function ldl_ltsolve!(n, x, Lp, Li, Lx)
  @inbounds for j = n:-1:1
    xj = x[j]
    @inbounds for p = Lp[j] : (Lp[j+1] - 1)
      xj -= Lx[p] * x[Li[p]]
    end
    x[j] = xj
  end
  return x
end

function ldl_solve(n, b, Lp, Li, Lx, D, P)
  y = b[P]
  ldl_lsolve!(n, y, Lp, Li, Lx)
  ldl_dsolve!(n, y, D)
  ldl_ltsolve!(n, y, Lp, Li, Lx)
  x = similar(b)
  x[P] = y
  return x
end

# a simplistic type for LDL' factorizations so we can do \
mutable struct LDLFactorization{T<:Real,Ti<:Integer}
  L::SparseMatrixCSC{T,Ti}
  D::Vector{T}
  P::Vector{Ti}
end

# use AMD permutation by default
ldl(A::Array{T,2}, args...) where T<:Real = ldl(sparse(A), args...)
ldl(A::SparseMatrixCSC{T,Ti}) where {T<:Real,Ti<:Integer} = ldl(A, amd(A))

# use ldl(A, collect(1:n)) to suppress permutation
function ldl(A::SparseMatrixCSC{T,Ti}, P::Vector{Int}) where {T<:Real,Ti<:Integer}
  n = size(A, 1)
  n == size(A, 2) || throw(DimensionMismatch("matrix must be square"))
  n == length(P) || throw(DimensionMismatch("permutation size mismatch"))

  # allocate space for symbolic analysis
  parent = Vector{Ti}(undef, n)
  Lnz = Vector{Ti}(undef, n)
  flag = Vector{Ti}(undef, n)
  pinv = Vector{Ti}(undef, n)
  Lp = Vector{Ti}(undef, n+1)

  # perform symbolic analysis
  ldl_symbolic!(n, A.colptr, A.rowval, Lp, parent, Lnz, flag, P, pinv)

  # allocate space for numerical factorization
  Li = Vector{Ti}(undef, Lp[n] - 1)
  Lx = Vector{T}(undef, Lp[n] - 1)
  Y = Vector{T}(undef, n)
  D = Vector{T}(undef, n)
  pattern = Vector{Ti}(undef, n)

  # perform numerical factorization
  ldl_numeric!(n, A.colptr, A.rowval, A.nzval, Lp, parent, Lnz, Li, Lx, D, Y, pattern, flag, P, pinv)

  return LDLFactorization(SparseMatrixCSC{T,Ti}(n, n, Lp, Li, Lx), D, P)
end

import Base.(\)
(\)(LDL::LDLFactorization{T,Ti}, b::AbstractVector{T}) where {T<:Real,Ti<:Integer} =
  ldl_solve(LDL.L.n, b, LDL.L.colptr, LDL.L.rowval, LDL.L.nzval, LDL.D, LDL.P)

end  # module
