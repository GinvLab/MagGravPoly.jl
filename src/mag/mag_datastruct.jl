
##########################################

"""
$(TYPEDEF)

Structure containing the components of a magnetization vector, 
  i.e., module, inclination and declination angles.

# Fields 

$(TYPEDFIELDS)

"""
Base.@kwdef struct MagnetizVector
    "modulus"
    mod::AbstractVector{<:Real}
    "inclination in degrees"
    Ideg::AbstractVector{<:Real}
    "declination in degrees"
    Ddeg::AbstractVector{<:Real}
end

##########################################

"""
$(TYPEDEF)

Structure containing a set of polygonal bodies (described by their segments and all vertices) along with their magnetizations (Induced + Remanent).
To create an instance, input an array of vectors of indices 
  (of vertices) for each body and the array of all the vertices.

# Fields 

$(TYPEDFIELDS)

"""
struct MagPolygBodies2D
    "structure defining the geometry of the bodies"
    geom::PolygBodies2D
    "vector of induced magnetizations"
    Jind::MagnetizVector
    "vector of remnant magnetizations"
    Jrem::MagnetizVector
    "[y1 y2] = lateral extension of the polygonal body. Tipycally y1 is negative since observations are at y=0."
    ylatext::Union{Vector{<:Real},Nothing}
    
    function MagPolygBodies2D(bodyindices::Vector{<:Vector{<:Integer}},allvert::AbstractArray{<:Real,2},
                              Jind::MagnetizVector,Jrem::MagnetizVector ;
                              ylatext::Union{Real,Vector{<:Real},Nothing} )
                              #y1::Union{Real,Nothing}, y2::Union{Real,Nothing}=nothing )
        ## define the geometry of bodies
        geom = PolygBodies2D(bodyindices,allvert)
        ## magnetic properties
        N=length(bodyindices)
        @assert (length(Jind.mod)==N && length(Jrem.mod)==N)
        @assert (length(Jind.Ideg)==N && length(Jrem.Ideg)==N)
        @assert (length(Jind.Ddeg)==N && length(Jrem.Ddeg)==N)
        # # induced magnetization
        # Jind = MagnetizVector(mod=modind,Ideg=Idegind,Ddeg=Ddegind)
        # # remnant magnetization
        # Jrem = MagnetizVector(mod=modrem,Ideg=Idegrem,Ddeg=Ddegrem)

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
            error("MagPolygBodies2D(): the following options are the only valid to specify:
                            - ylatext=`nothing` for the pure 2D formulation;
                            - specify only a Real number for ylatext (the code then assumes -y1=y2) for the 2.5D formulation (symmetric case);
                            - specify both y1 and y2 as ylatext=[y1, y2] for the 2.75D formulation (asymmetric case).
                            Aborting!")
        end

        return new(geom,Jind,Jrem,outylatext)
    end
end

########################################################################
