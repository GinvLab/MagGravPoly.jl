

using Documenter, MagGravPoly


makedocs(repo="https://github.com/GinvLab/MagGravPoly.jl",
         sitename="MagGravPoly.jl",
         modules = [MagGravPoly],
         authors = "Alessandro Ghirotto, Andrea Zunino",
         format = Documenter.HTML(prettyurls=get(ENV,"CI",nothing)=="true"),
         pages = [
             "Home - MagGravPoly" => "index.md",
             "GeoPoly" => "geopoly.md"
         ],
         warnonly = [:missing_docs, :cross_references]
         )

deploydocs(
    repo="github.com/GinvLab/MagGravPoly.jl.git",
    devbranch = "main",
    deploy_config = Documenter.GitHubActions(),
    branch = "gl-pages"
)



###########################################################
