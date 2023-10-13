

using Documenter, MagGrav2Dpoly


makedocs(repo="https://gitlab.com/JuliaGeoph/MagGrav2Dpoly.jl/blob/{commit}{path}#{line}",
         sitename="MagGrav2Dpoly.jl",
         modules = [MagGrav2Dpoly],
         authors = "Andrea Zunino, Alessandro Ghirotto",
         format = Documenter.HTML(prettyurls=get(ENV,"CI",nothing)=="true"),
         pages = [
             "Home" => "index.md",
         ],
         warnonly = [:missing_docs, :cross_references]
         )

deploydocs(
    repo="gitlab.com/JuliaGeoph/MagGrav2Dpoly.jl.git",
    devbranch = "main",
    deploy_config = Documenter.GitLab(),
    branch = "gl-pages"
)



###########################################################
