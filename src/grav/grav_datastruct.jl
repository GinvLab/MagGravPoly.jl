
##########################################

"""
$(TYPEDEF)

Structure containing a set of polygonal bodies (described by their segments and all vertices) along with their densities.
To create an instance, input an array of vectors of indices (of vertices) for each body and the array of all the vertices.

# Fields 

$(TYPEDFIELDS)

"""
struct GravPolygBodies2D
    "structure defining the geometry of the bodies"
    geom::PolygBodies2D
    "vector of densities"
    rho::Vector{<:Real}
    "[y1 y2] = lateral extension of the polygonal body. Tipycally y1 is negative since observations are at y=0."
    ylatext::Union{Vector{<:Real},Nothing}

    function GravPolygBodies2D(bodyindices::Vector{<:Vector{<:Integer}},allvert::AbstractArray{<:Real,2},
                               rho::Vector{<:Real};
                               ylatext::Union{Real,Vector{<:Real},Nothing} )
        
        ## define the geometry of bodies
        geom = PolygBodies2D(bodyindices,allvert)
        ## density vector dim check
        N=length(bodyindices)
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
            error("GravPolygBodies2D(): the following options are the only valid to specify:
                                - ylatext=`nothing` for the pure 2D formulation;
                                - specify only a Real number for ylatext (the code then assumes -y1=y2) for the 2.5D formulation (symmetric case);
                                - specify both y1 and y2 as ylatext=[y1, y2] for the 2.75D formulation (asymmetric case).
                                Aborting!")
        end
        
        return new(geom,rho,outylatext)
    end
end

##########################################


















