##############################
"""
GeoPolygons

A module developed for handling polygonal shapes in the framework of 2D potential field forward and inverse modeling.

# Exports

$(EXPORTS)
"""
module GeoPolygons

using DocStringExtensions
using LibGEOS
using Random
using Statistics

export PolygBodies2D,BodySegments2D,TopoEdges
export checktopo,checkpoly,selfintersectpoly
export checkanticlockwiseorder,checkbodyindices
export calcareapoly,calcareamanypoly


include("polydatastruct.jl")
include("interseg.jl")
include("polygtests.jl")
include("topogtests.jl")
include("checkutils.jl")
include("utils.jl")

end # module
