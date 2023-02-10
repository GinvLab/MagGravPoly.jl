

using Documenter, JointMagGrav2Dpoly


makedocs(sitename="MagGrav2Dpoly.jl",
         modules = [MagGrav2Dpoly],
         authors = "Andrea Zunino, Alessandro Ghirotto",
         format = Documenter.HTML(prettyurls=get(ENV,"CI",nothing)=="true"),
         pages = [
             "Home" => "index.md",
         ]
         )

deploydocs(
    repo="gitlab.com/JuliaGeoph/MagGrav2Dpoly.jl.git",
    devbranch = "main",
    deploy_config = Documenter.GitLab(),
    branch = "gl-pages"
)



###########################################################
