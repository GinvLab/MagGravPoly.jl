#################################################################

l12(n,ncoo,nbo) = (ncoo+1+(n-1)*nbo,ncoo+n*nbo)

#################################################################
"""
$(TYPEDSIGNATURES)

Function to create a vector of model parameters the user would like to invert from a `MagPolygBodies2D` polygons structure. For this purpose is required as input even a `Mag2DPolyMisf` misfit structure.  
"""
function magstruct2vec(whichpar::Symbol,magpbod::MagPolygBodies2D)
    
    if whichpar==:all

        nbo = magpbod.geom.nbo
        nvert = size(magpbod.geom.allvert,1)
        ncoo = Integer(nvert*2)
        modpar = zeros(eltype(magpbod.geom.allvert),ncoo+(6*nbo))

        # set vertices
        modpar[1:ncoo] .= magpbod.geom.allvert[:]

        #l=ncoo+1
        # set induced magnetization
        l1,l2 = l12(1,ncoo,nbo)
        modpar[l1:l2] .= magpbod.Jind.mod[:]
        l1,l2 = l12(2,ncoo,nbo)
        modpar[l1:l2] .= magpbod.Jind.Ideg[:]
        l1,l2 = l12(3,ncoo,nbo)
        modpar[l1:l2] .= magpbod.Jind.Ddeg[:]

        # set remnant magnetization
        l1,l2 = l12(4,ncoo,nbo)
        modpar[l1:l2] .= magpbod.Jrem.mod[:]
        l1,l2 = l12(5,ncoo,nbo)
        modpar[l1:l2] .= magpbod.Jrem.Ideg[:]
        l1,l2 = l12(6,ncoo,nbo)
        modpar[l1:l2] .= magpbod.Jrem.Ddeg[:]

    elseif whichpar==:vertices

        nvert = size(magpbod.geom.allvert,1)
        ncoo = Integer(nvert*2)
        modpar = zeros(eltype(magpbod.geom.allvert),ncoo)

        #set only vertices
        modpar[1:ncoo] .= magpbod.geom.allvert[:]
        
    elseif whichpar==:magnetization

        nbo = magpbod.geom.nbo
        modpar = zeros(eltype(magpbod.geom.allvert),nbo*6)

        # set induced magnetization
        modpar[1:nbo] .= magpbod.Jind.mod[:]
        modpar[nbo+1:2*nbo] .= magpbod.Jind.Ideg[:]
        modpar[(2*nbo)+1:3*nbo] .= magpbod.Jind.Ddeg[:]

        # set remnant magnetization
        modpar[(3*nbo)+1:4*nbo] .= magpbod.Jrem.mod[:]
        modpar[(4*nbo)+1:5*nbo] .= magpbod.Jrem.Ideg[:]
        modpar[(5*nbo)+1:6*nbo] .= magpbod.Jrem.Ddeg[:]
        
    end      
        
    return modpar
end

##########################################
#=
function copyintomagstruct!(magpbod::MagPolygBodies2D,modpar::AbstractVector{<:Real})

    nbo = magpbod.geom.nbo
    ncoo = length(modpar)-6*nbo
    nvert = div(ncoo,2)

    # set vertices
    allvert = reshape(modpar[1:ncoo],nvert,2)
    @show size(modpar)
    @show ncoo,nvert
    @show size(allvert)
    @show size(magpbod.geom.allvert)
    magpbod.geom.allvert .= allvert

    # set magnetization
    l12(n,ncoo,nbo) = (ncoo+1+(n-1)*nbo,ncoo+1+n*nbo-1)

    l1m,l2m = l12(1,ncoo,nbo)
    l1I,l2I = l12(2,ncoo,nbo)
    l1D,l2D = l12(3,ncoo,nbo)
    magpbod.Jind.mod  .= modpar[l1m:l2m]
    magpbod.Jind.Ideg .= modpar[l1I:l2I]
    magpbod.Jind.Ddeg .= modpar[l1D:l2D]
    
    l1m,l2m = l12(4,ncoo,nbo)
    l1I,l2I = l12(5,ncoo,nbo)
    l1D,l2D = l12(6,ncoo,nbo)
    magpbod.Jrem.mod  .= modpar[l1m:l2m]
    magpbod.Jrem.Ideg .= modpar[l1I:l2I]
    magpbod.Jrem.Ddeg .= modpar[l1D:l2D]
    
    return magpbod
end
=#
##########################################
"""
$(TYPEDSIGNATURES)

Function to reconstruct a `MagPolygBodies2D` polygons structure from i) a list of `bodyindices` and ii) a vector of model parameters the user would like to invert. For this purpose is required as input even a `Mag2DPolyMisf` misfit structure.  
"""
function vecmodpar2magstruct(magmisf::Mag2DPolyMisf,curbodyindices::Vector{<:Vector{<:Integer}},modpar::Vector{<:Real})
    ## separate "curbodyindices" is necessary because the bodyindices might not be in sync
    ##  with the iterations of HMC, so they must be taken from the .h5 file.

    if magmisf.whichpar==:all

        nbo = length(magmisf.bodyindices)
        ncoo = length(modpar)-(6*nbo)
        nvert = div(ncoo,2)
        
        # set vertices
        allvert = reshape(modpar[1:ncoo],nvert,2)

        # set induced magnetization
        l1m,l2m = l12(1,ncoo,nbo)
        l1I,l2I = l12(2,ncoo,nbo)
        l1D,l2D = l12(3,ncoo,nbo)
        Jind = MagnetizVector(mod=modpar[l1m:l2m], Ideg=modpar[l1I:l2I],
                              Ddeg=modpar[l1D:l2D])

        # set remnant magnetization
        l1m,l2m = l12(4,ncoo,nbo)
        l1I,l2I = l12(5,ncoo,nbo)
        l1D,l2D = l12(6,ncoo,nbo)
        Jrem = MagnetizVector(mod=modpar[l1m:l2m], Ideg=modpar[l1I:l2I],
                              Ddeg=modpar[l1D:l2D])
        
    elseif magmisf.whichpar==:vertices

        ncoo = length(modpar) 
        nvert = div(ncoo,2)
        
        # set vertices
        allvert = reshape(modpar[1:ncoo],nvert,2)

        # set magnetization
        Jind = deepcopy(magmisf.Jind)
        Jrem = deepcopy(magmisf.Jrem)
        
    elseif magmisf.whichpar==:magnetization

        nbo = length(magmisf.bodyindices)
        
        # set vertices
        allvert = magmisf.allvert

        # set induced magnetization
        Jind_mod = modpar[1:nbo]
        Jind_Ideg = modpar[nbo+1:2*nbo]
        Jind_Ddeg = modpar[(2*nbo)+1:3*nbo]
        Jind = MagnetizVector(mod=Jind_mod[1:end], Ideg=Jind_Ideg[1:end],
                              Ddeg=Jind_Ddeg[1:end])

        # set remnant magnetization
        Jrem_mod = modpar[(3*nbo)+1:4*nbo]
        Jrem_Ideg = modpar[(4*nbo)+1:5*nbo]
        Jrem_Ddeg = modpar[(5*nbo)+1:6*nbo]        
        Jrem = MagnetizVector(mod=Jrem_mod[1:end], Ideg=Jrem_Ideg[1:end],
                              Ddeg=Jrem_Ddeg[1:end])
        
    end

    mbody = MagPolygBodies2D(curbodyindices,allvert,Jind,Jrem,
                             ylatext=magmisf.ylatext)

    return mbody
end

#############################################
##########################################

"""
$(TYPEDSIGNATURES)

Function to extract a list of polygon vertices from i) a `Mag2DPolyMisf` misfit structure and ii) a vector of model parameters the user would like to invert.  
"""
function vecmodpar2magvertices(magmisf::Mag2DPolyMisf,modpar::Vector{<:Real})

    whichparmag = magmisf.whichpar

    nbo = length(magmisf.bodyindices)

    if whichparmag==:all 
        
        ncoo = length(modpar)-(6*nbo)
        nvert = div(ncoo,2)

        # set vertices
        allvert = reshape(modpar[1:ncoo],nvert,2)


    elseif whichparmag==:vertices 
        
        ncoo = length(modpar)##-nbo
        nvert = div(ncoo,2)

        # set vertices
        allvert = reshape(modpar[1:ncoo],nvert,2)
        
    elseif whichparmag==:magnetization 

        # set vertices
        allvert = copy(magmisf.allvert)

    else
        error("vecmodpar2magvertices(): Wrong argument 'magmisf.whichpar'. Aborting.")
    end
    
    return allvert
end

############################################################################
