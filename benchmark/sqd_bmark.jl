using BenchmarkTools
using MatrixMarket

using LinearAlgebra
using Printf
using SparseArrays

using LDLFactorizations
using QDLDL

# download from https://github.com/optimizers/sqd-collection
const sqd_path = "/Users/dpo/dev/datasets/sqd-collection"
subdirs = readdir(sqd_path)
const formulations = ("2x2", "3x3")
const iters = (0, 5, 10)

const SUITE = BenchmarkGroup()

for subdir ∈ subdirs
  subdir == ".git" && continue
  isdir(joinpath(sqd_path, subdir)) || continue  # ignore regular files
  for formulation ∈ formulations
    for iter ∈ iters
      iterpath = joinpath(sqd_path, subdir, formulation, "iter_$(iter)")
      isdir(iterpath) || continue
      A = MatrixMarket.mmread(joinpath(iterpath, "K_$(iter).mtx"))
      name = "$(subdir)_$(formulation)_$(iter)"
      SUITE[name] = BenchmarkGroup()
      SUITE[name]["LDL"] = @benchmarkable ldl($A)
      SUITE[name]["QDLDL"] = @benchmarkable qdldl($A)
    end
  end
end
