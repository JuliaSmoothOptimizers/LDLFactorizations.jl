using DataFrames
using GitHub
using JLD2
using JSON
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

# NB: benchmarkpkg will run benchmarks/benchmarks.jl by default
commit = benchmarkpkg("LDLFactorizations")  # current state of repository
master = benchmarkpkg("LDLFactorizations", "master")
judgement = judge(commit, master)

# TODO: automate the following
commit_stats = bmark_results_to_dataframes(commit)
master_stats = bmark_results_to_dataframes(master)
judgement_stats = judgement_results_to_dataframes(judgement)

# extract stats for each benchmark to plot profiles
fact_stats = Dict{Symbol,DataFrame}(:commit => commit_stats[:fact],
                                    :master => master_stats[:fact])
solve1_stats = Dict{Symbol,DataFrame}(:commit => commit_stats[:solve1],
                                      :master => master_stats[:solve1])
solve5_stats = Dict{Symbol,DataFrame}(:commit => commit_stats[:solve5],
                                      :master => master_stats[:solve5])
save_stats(fact_stats, "ldl_commit_vs_master_fact.jld2")
save_stats(solve1_stats, "ldl_commit_vs_master_solve1.jld2")
save_stats(solve5_stats, "ldl_commit_vs_master_solve5.jld2")
save_stats(judgement_stats, "ldl_commit_vs_master_judgement.jld2")

export_markdown("judgement.md", judgement)
export_markdown("master.md", master)
export_markdown("commit.md", commit)

function profile_solvers_from_pkgbmark(stats::Dict{Symbol,DataFrame})
  # guard against zero gctimes
  costs = [df -> df[!, :time], df -> df[!, :memory], df -> df[!, :gctime] .+ 1, df -> df[!, :allocations]]
  profile_solvers(stats, costs, ["time", "memory", "gctime+1", "allocations"])
end

# save and read profiles in SVG format
p_fact = profile_solvers_from_pkgbmark(fact_stats)
savefig("profiles_commit_vs_master_fact.svg")
svgfile_fact = open("profiles_commit_vs_master_fact.svg", "r") do fd
  readlines(fd)
end

p_solve1 = profile_solvers_from_pkgbmark(solve1_stats)
savefig("profiles_commit_vs_master_solve1.svg")
svgfile_solve1 = open("profiles_commit_vs_master_solve1.svg", "r") do fd
  readlines(fd)
end

p_solve5 = profile_solvers_from_pkgbmark(solve5_stats)
savefig("profiles_commit_vs_master_solve5.svg")
svgfile_solve5 = open("profiles_commit_vs_master_solve5.svg", "r") do fd
  readlines(fd)
end

gist_json = JSON.parse("""
    {
      "description": "LDLFactorizations repository benchmark",
      "public": true,
      "files": {
        "1_fact.svg": {
          "content": "$(escape_string(join(svgfile_fact)))"
        },
        "2_solve1.svg": {
          "content": "$(escape_string(join(svgfile_solve1)))"
        },
        "3_solve5.svg": {
          "content": "$(escape_string(join(svgfile_solve5)))"
        },
        "4_judgement.md": {
          "content": "$(escape_string(sprint(export_markdown, judgement)))"
        },
        "5_commit.md": {
          "content": "$(escape_string(sprint(export_markdown, commit)))"
        },
        "6_master.md": {
          "content": "$(escape_string(sprint(export_markdown, master)))"
        }
      }
    }""")

# Need to add GITHUB_AUTH to your .bashrc
myauth = GitHub.authenticate(ENV["GITHUB_AUTH"])
posted_gist = create_gist(params = gist_json, auth = myauth)
println(posted_gist.html_url)
