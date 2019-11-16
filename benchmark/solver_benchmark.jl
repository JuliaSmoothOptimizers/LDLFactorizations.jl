using DataFrames
using JLD2
using PkgBenchmark
using Plots

using SolverBenchmark

pyplot()
using PlotThemes
theme(:wong)

common_plot_args = Dict{Symbol,Any}(
  :linewidth => 2,
  :alpha => .75,
  :titlefontsize => 8,
  :legendfontsize => 8,
  :xtickfontsize => 6,
  :ytickfontsize => 6,
  :guidefontsize => 8,
)
Plots.default(; common_plot_args...)

function save_stats(stats::Dict{Symbol,DataFrame}, filename::AbstractString; force::Bool=false, key::String="stats")
  isfile(filename) && !force && error("$filename already exists; use `force=true` to overwrite")
  jldopen(filename, "w") do file
    file[key] = stats
  end
end

# perform benchmarks
results = PkgBenchmark.benchmarkpkg("LDLFactorizations", script="benchmark/sqd_bmark.jl")
stats = bmark_results_to_dataframes(results)
save_stats(stats, "ldl_vs_qdldl.jld2")

# process benchmark results and post gist
p = profile_solvers(results)
posted_gist = to_gist(results, p)
println(posted_gist.html_url)
