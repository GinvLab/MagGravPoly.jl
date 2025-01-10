

"""
MagGravPoly

"""
module MagGravPoly

#export GeoPoly
export MG2D

# include but do not export GeoPoly
include("GeoPoly/GeoPoly.jl")

# MagGrav2Dpoly
include("MG2D/MG2D.jl")


end # module MagGravPoly
