
"""
MagGrav2Dpoly

A module to perform magnetic and gravity anomaly calculations for 2D polygonal bodies.

# Exports

$(EXPORTS)
"""
module MagGrav2Dpoly


using GeoPolygons
using LinearAlgebra
using ForwardDiff
using ReverseDiff
using DocStringExtensions

#using Statistics
using PrettyTables

#####################################
## Grav-related stuff

export GravPolygBodies2D
export tgravpolybodies2D,tgravpolybodies2Dgen

export Grav2DPolyMisf 
export calcmisfgrav,precalcADstuffgrav,calc∇misfgrav

export vecmodpar2gravstruct,gravstruct2vec
export calcMgrav,printgravmodparinfo

include("grav/grav_datastruct.jl")
include("grav/grav_derivatives.jl")
include("grav/grav_struct2vec.jl")
include("grav/grav_2dpolybodies.jl")
include("grav/grav_utils.jl")

####################################

#####################################
## Mag-related stuff

export MagPolygBodies2D
export tmagpolybodies2D,tmagpolybodies2Dgen
export MagnetizVector

export Mag2DPolyMisf 
export calcmisfmag,precalcADstuffmag,calc∇misfmag

export vecmodpar2magstruct,magstruct2vec
export calcMmag,printmagmodparinfo

include("mag/mag_datastruct.jl")
include("mag/mag_derivatives.jl")
include("mag/mag_struct2vec.jl")
include("mag/mag_2dpolybodies.jl")
include("mag/mag_utils.jl")

####################################


#####################################
## Joint Mag Grav-related stuff

export JointPolygBodies2D
export tjointpolybodies2D,tjointpolybodies2Dgen

export Joint2DPolyMisf 
export calcmisfjointmaggrav,precalcADstuffjointmaggrav,calc∇misfjointmaggrav

export vecmodpar2jointstruct,jointstruct2vec
export calcMjoint,printjointmodparinfo

include("joint/joint_datastruct.jl")
include("joint/joint_derivatives.jl")
include("joint/joint_struct2vec.jl")
include("joint/joint_2dpolybodies.jl")
include("joint/joint_utils.jl")

####################################


####################################
## HMC stuff

# grav
include("grav/HMCgrav2Dpoly.jl")
using .HMCGrav2Dpoly
export Grav2DpolyProb,gravprob2vec

# mag
include("mag/HMCmag2Dpoly.jl")
using .HMCMag2Dpoly
export Mag2DpolyProb,magprob2vec

# joint
include("joint/HMCJointMagGrav2Dpoly.jl")
using .HMCJointMagGrav2Dpoly
export Joint2DpolyProb,jointprob2vec


####################################




end # module
