###################################################################
"""
$(TYPEDSIGNATURES)

Function to create a vector of model parameters the user would like to invert from a `GravPolygBodies2D` polygons structure. 
For this purpose is required as input even a `Grav2DPolyMisf` misfit structure.  
"""
function gravstruct2vec(whichpar::Symbol,gravpbod::GravPolygBodies2D)

    if whichpar==:all
        
        nbo = gravpbod.geom.nbo
        nvert = size(gravpbod.geom.allvert,1)
        ncoo = Integer(nvert*2)
        modpar = zeros(eltype(gravpbod.geom.allvert),ncoo+nbo)

        #set vertices
        modpar[1:ncoo] .= gravpbod.geom.allvert[:]
        # set density
        modpar[ncoo+1:end] .= gravpbod.rho[:]

    elseif whichpar==:vertices

        nvert = size(gravpbod.geom.allvert,1)
        ncoo = Integer(nvert*2)
        modpar = zeros(eltype(gravpbod.geom.allvert),ncoo)

        # set only vertices
        modpar[1:ncoo] .= gravpbod.geom.allvert[:]

    elseif whichpar==:density

        nbo = gravpbod.geom.nbo
        modpar = zeros(eltype(gravpbod.geom.allvert),nbo)

        # set only density
        modpar[1:nbo] .= gravpbod.rho[:]
        
    end
            
    return modpar
end

##########################################
#=
function copyintogravstruct!(gravpbod::GravPolygBodies2D,modpar::AbstractVector{<:Real})

    nbo = gravpbod.geom.nbo
    ncoo = length(modpar)-nbo
    nvert = div(ncoo,2)

    # set vertices
    allvert = reshape(modpar[1:ncoo],nvert,2)
    @show size(modpar)
    @show ncoo,nvert
    @show size(allvert)
    @show size(gravpbod.geom.allvert)
    gravpbod.geom.allvert .= allvert

    # set densities
    l12(n,ncoo,nbo) = (ncoo+1+(n-1)*nbo,ncoo+nbo+(n-1)*nbo)

    l1d,l2d = l12(1,ncoo,nbo)  
    gravpbod.rho .= modpar[l1d:l2d]
    
    return gravpbod
end
=#
##########################################
"""
$(TYPEDSIGNATURES)

Function to reconstruct a `GravPolygBodies2D` polygons structure from the vector of model parameters the user would like to invert. 
for this purpose is required as input even a `Grav2DPolyMisf` misfit structure.  
"""
function vecmodpar2gravstruct(gravmisf::Grav2DPolyMisf,modpar::Vector{<:Real})

    if gravmisf.whichpar==:all

        nbo = length(gravmisf.bodyindices)
        ncoo = length(modpar)-nbo
        nvert = div(ncoo,2)
        
        # set vertices
        allvert = reshape(modpar[1:ncoo],nvert,2)
        # set density
        rho = modpar[ncoo+1:end]

    elseif gravmisf.whichpar==:vertices

        ncoo = length(modpar) 
        nvert = div(ncoo,2)

        # set vertices
        allvert = reshape(modpar[1:ncoo],nvert,2)     
        # set density
        rho = gravmisf.rho

    elseif gravmisf.whichpar==:density

        # set vertices
        allvert = gravmisf.allvert
        # set density
        rho = modpar[1:end]
                
    end       
    
    mbody = GravPolygBodies2D(gravmisf.bodyindices,allvert,rho,
                              ylatext=gravmisf.ylatext)

    return mbody
end

#############################################
