

"""
MagGravPoly

"""
module MagGravPoly

#export GeoPolygons
export MG2D

# include but do not export GeoPolygons
include("GeoPolygons/GeoPolygons.jl")

# MagGrav2Dpoly
include("MG2D/MG2D.jl")


end # module MagGravPoly
