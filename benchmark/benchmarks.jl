# benchmark file to compare two commits of LDLFactorizations

using BenchmarkTools
using MatrixMarket

using DelimitedFiles
using LinearAlgebra
using Pkg.Artifacts
using Printf
using SparseArrays

using LDLFactorizations

# obtain path to SQD collection
const artifact_toml = joinpath(@__DIR__, "Artifacts.toml")
const sqd_hash = artifact_hash("sqdcollection", artifact_toml)
@assert artifact_exists(sqd_hash)
const sqd_path = joinpath(artifact_path(sqd_hash), "sqd-collection-0.1")
subdirs = readdir(sqd_path)
const formulations = ("2x2", "3x3")
const iters = (0, 5, 10)

const SUITE = BenchmarkGroup()
SUITE["fact"] = BenchmarkGroup()
SUITE["1solve"] = BenchmarkGroup()
SUITE["5solve"] = BenchmarkGroup()

for subdir ∈ subdirs
  subdir == ".git" && continue
  isdir(joinpath(sqd_path, subdir)) || continue  # ignore regular files
  for formulation ∈ formulations
    for iter ∈ iters
      iterpath = joinpath(sqd_path, subdir, formulation, "iter_$(iter)")
      isdir(iterpath) || continue
      A = MatrixMarket.mmread(joinpath(iterpath, "K_$(iter).mtx"))
      b = readdlm(joinpath(iterpath, "rhs_$(iter).rhs"))[:, 1]
      B = [b b b b b]
      name = "$(subdir)_$(formulation)_$(iter)"
      SUITE["fact"][name] = @benchmarkable ldl($A)
      LDL = ldl(A)
      x = similar(b)
      SUITE["1solve"][name] = @benchmarkable ldiv!($x, $LDL, $b)
      X = similar(B)
      SUITE["5solve"][name] = @benchmarkable ldiv!($X, $LDL, $B)
    end
  end
end
