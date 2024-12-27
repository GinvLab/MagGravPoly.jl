
#############################################
"""
$(TYPEDEF)

Structure containing a set of polygonal bodies (described by their segments and all vertices) along with their magnetizations (Induced + Remanent) and densities.
To create an instance, input an array of vectors of indices 
  (of vertices) for each body and the array of all the vertices.

# Fields 

$(TYPEDFIELDS)

"""
Base.@kwdef struct JointPolygBodies2D
    "structure defining the geometry of the bodies"
    geom::PolygBodies2D
    "vector of induced magnetizations"
    Jind::MagnetizVector
    "vector of remnant magnetizations"
    Jrem::MagnetizVector
    "vector of densities"
    rho::Vector{<:Real}
    "[y1 y2] = lateral extension of the polygonal body. Tipycally y1 is negative since observations are at y=0."
    ylatext::Union{Vector{<:Real},Nothing}
    
    function JointPolygBodies2D(bodyindices::Vector{<:Vector{<:Integer}},allvert::AbstractArray{<:Real,2},
                                Jind::MagnetizVector,Jrem::MagnetizVector,rho::Vector{<:Real};
                                ylatext::Union{Real,Vector{<:Real},Nothing} )
        
        ## define the geometry of bodies
        geom = PolygBodies2D(bodyindices,allvert)
        ## magnetic properties
        N=length(bodyindices)
        @assert (length(Jind.mod)==N && length(Jrem.mod)==N)
        @assert (length(Jind.Ideg)==N && length(Jrem.Ideg)==N)
        @assert (length(Jind.Ddeg)==N && length(Jrem.Ddeg)==N)
        ## density vector dim check
        @assert length(rho)==N
        
        ## Optional 2.5/2.75D dimensions
        if ylatext==nothing 
            outylatext = nothing

        elseif typeof(ylatext)<:Real 
            @assert ylatext > 0.0
            outylatext = [-ylatext, ylatext]

        elseif typeof(ylatext)<:Vector{<:Real} 
            @assert ylatext[1] < ylatext[2] 
            outylatext = copy(ylatext)

        else
            error("JointPolygBodies2D(): the following options are the only valid to specify:
                                - ylatext=`nothing` for the pure 2D formulation;
                                - specify only a Real number for ylatext (the code then assumes -y1=y2) for the 2.5D formulation (symmetric case);
                                - specify both y1 and y2 as ylatext=[y1, y2] for the 2.75D formulation (asymmetric case).
                                Aborting!")
        end
        
        return new(geom,Jind,Jrem,rho,outylatext)
    end
end

#################################################################
# """
#  Function to perform merging between MagPolygBodies2D and GravPolygBodies2D structures in a JointPolygBodies2D one.
# """
# function MergePolygBodies2D(magpbod::MagPolygBodies2D,gravpbod::GravPolygBodies2D)

#     # check bodies are the same between mag and grav structures
#     @assert length(magpbod.geom.allvert) == length(gravpbod.geom.allvert)
#     @assert all(magpbod.geom.allvert .== gravpbod.geom.allvert)
#     @assert all(magpbod.geom.bodyindices .== gravpbod.geom.bodyindices)
#     @assert length(magpbod.Jind.mod) == length(gravpbod.rho)
    
#     ## create a structure containing geometries of the bodies and their mag and grav properties
#     jointpbod = JointPolygBodies2D(magpbod.geom.bodyindices,magpbod.geom.allvert,magpbod.Jind,magpbod.Jrem,gravpbod.rho)
 
#     return jointpbod
# end

#################################################################
