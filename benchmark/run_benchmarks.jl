using GitHub, JSON, PkgBenchmark

# NB: benchmarkpkg will run benchmarks/benchmarks.jl by default
commit = benchmarkpkg("LDLFactorizations")  # current state of repository
master = benchmarkpkg("LDLFactorizations", "master")
judgement = judge(commit, master)
export_markdown("judgement.md", judgement)
export_markdown("master.md", master)
export_markdown("commit.md", commit)
open("profiles.svg", "r") do fd
  svgfile = readlines(fd)
end

gist_json = JSON.parse("""
    {
        "description": "LDLFactorizations repository benchmark",
        "public": true,
        "files": {
            "judgement.md": {
                "content": "$(escape_string(sprint(export_markdown, judgement)))"
            },
            "master.md": {
                "content": "$(escape_string(sprint(export_markdown, master)))"
            },
            "commit.md": {
                "content": "$(escape_string(sprint(export_markdown, commit)))"
            },
            "profile.svg": {
                "content": "$(escape_string(svgfile))"
            }
        }
    }""")

# Need to add GITHUB_AUTH to your .bashrc
myauth = GitHub.authenticate(ENV["GITHUB_AUTH"])
posted_gist = create_gist(params = gist_json, auth = myauth)
println(posted_gist.html_url)
