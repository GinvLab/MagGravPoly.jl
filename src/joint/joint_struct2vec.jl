#######################################################################

## already defined in the mag files
## l12(n,ncoo,nbo) = (ncoo+1+(n-1)*nbo,ncoo+n*nbo)

#######################################################################
"""
$(TYPEDSIGNATURES)

Function to create a vector of model parameters proving as input i) the type of parameters the user would like to invert for (gravity + magnetic problems) and ii) a `JointPolygBodies2D` polygons structure.  
"""
function jointstruct2vec(mag_whichpar::Symbol,grav_whichpar::Symbol,jointpbod::JointPolygBodies2D)

    whichparmag = mag_whichpar
    whichpargrav = grav_whichpar
    
    if whichparmag==:all && whichpargrav==:all
        
        nbo = jointpbod.geom.nbo
        nvert = size(jointpbod.geom.allvert,1)
        ncoo = Integer(nvert*2)
        modpar = zeros(eltype(jointpbod.geom.allvert),ncoo+(7*nbo))

        # set vertices
        modpar[1:2*nvert] .= jointpbod.geom.allvert[:]

        # set induced magnetization
        #l=ncoo+1       
        l1,l2 = l12(1,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.Jind.mod[:]
        l1,l2 = l12(2,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.Jind.Ideg[:]
        l1,l2 = l12(3,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.Jind.Ddeg[:]

        # set remnant magnetization
        l1,l2 = l12(4,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.Jrem.mod[:]
        l1,l2 = l12(5,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.Jrem.Ideg[:]
        l1,l2 = l12(6,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.Jrem.Ddeg[:]

        # set density
        l1,l2 = l12(7,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.rho[:]

    elseif whichparmag==:all && whichpargrav==:vertices

        nbo = jointpbod.geom.nbo
        nvert = size(jointpbod.geom.allvert,1)
        ncoo = Integer(nvert*2)
        modpar = zeros(eltype(jointpbod.geom.allvert),ncoo+(6*nbo))

        # set vertices
        modpar[1:2*nvert] .= jointpbod.geom.allvert[:]

        # set induced magnetization
        #l=ncoo+1       
        l1,l2 = l12(1,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.Jind.mod[:]
        l1,l2 = l12(2,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.Jind.Ideg[:]
        l1,l2 = l12(3,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.Jind.Ddeg[:]

        # set remnant magnetization
        l1,l2 = l12(4,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.Jrem.mod[:]
        l1,l2 = l12(5,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.Jrem.Ideg[:]
        l1,l2 = l12(6,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.Jrem.Ddeg[:]
        
    elseif whichparmag==:vertices && whichpargrav==:all
        
        nbo = jointpbod.geom.nbo
        nvert = size(jointpbod.geom.allvert,1)
        ncoo = Integer(nvert*2)
        modpar = zeros(eltype(jointpbod.geom.allvert),ncoo+nbo)

        # set vertices
        modpar[1:2*nvert] .= jointpbod.geom.allvert[:]

        # set density
        l1,l2 = l12(1,ncoo,nbo)
        modpar[l1:l2] .= jointpbod.rho[:]

    elseif whichparmag==:magnetization && whichpargrav==:density

        nbo = jointpbod.geom.nbo
        modpar = zeros(eltype(jointpbod.geom.allvert),nbo*7)

        # set induced magnetization
        modpar[1:nbo] .= jointpbod.Jind.mod[:]
        modpar[nbo+1:2*nbo] .= jointpbod.Jind.Ideg[:]
        modpar[(2*nbo)+1:3*nbo] .= jointpbod.Jind.Ddeg[:]

        # set remanant magnetization
        modpar[(3*nbo)+1:4*nbo] .= jointpbod.Jrem.mod[:]
        modpar[(4*nbo)+1:5*nbo] .= jointpbod.Jrem.Ideg[:]
        modpar[(5*nbo)+1:6*nbo] .= jointpbod.Jrem.Ddeg[:]

        # set density
        modpar[(6*nbo)+1:end] .= jointpbod.rho[:]
        
    elseif whichparmag==:vertices && whichpargrav==:vertices

        nvert = size(jointpbod.geom.allvert,1)
        ncoo = Integer(nvert*2)
        modpar = zeros(eltype(jointpbod.geom.allvert),ncoo)

        # set only vertices
        modpar[1:ncoo] .= jointpbod.geom.allvert[:]

    else
        error("jointstruct2vec(): Wrong argument 'mag_whichpar' and/or 'grav_whichpar': possible choices are: 
                          - mag_whichpar==:vertices && grav_whichpar==:vertices;
                          - mag_whichpar==:vertices && grav_whichpar==:all;
                          - mag_whichpar==:all && grav_whichpar==:all;
                          - mag_whichpar==:all && grav_whichpar==:vertices;
                          - mag_whichpar==:magnetization && grav_whichpar==:density.
                          Aborting!")
    
    end
    
    return modpar
end


##########################################
#=
function copyintojointstruct!(jointpbod::JointPolygBodies2D,modpar::AbstractVector{<:Real})

    nbo = jointpbod.geom.nbo
    ncoo = length(modpar)-7*nbo
    nvert = div(ncoo,2)

    # set vertices
    allvert = reshape(modpar[1:ncoo],nvert,2)
    @show size(modpar)
    @show ncoo,nvert
    @show size(allvert)
    @show size(jointpbod.geom.allvert)
    jointpbod.geom.allvert .= allvert

    # set all parameters
    l12(n,ncoo,nbo) = (ncoo+1+(n-1)*nbo,ncoo+nbo+(n-1)*nbo)

    l1m,l2m = l12(1,ncoo,nbo)
    l1I,l2I = l12(2,ncoo,nbo)
    l1D,l2D = l12(3,ncoo,nbo)
    jointpbod.Jind.mod  .= modpar[l1m:l2m]
    jointpbod.Jind.Ideg .= modpar[l1I:l2I]
    jointpbod.Jind.Ddeg .= modpar[l1D:l2D]
    
    l1m,l2m = l12(4,ncoo,nbo)
    l1I,l2I = l12(5,ncoo,nbo)
    l1D,l2D = l12(6,ncoo,nbo)
    jointpbod.Jrem.mod  .= modpar[l1m:l2m]
    jointpbod.Jrem.Ideg .= modpar[l1I:l2I]
    jointpbod.Jrem.Ddeg .= modpar[l1D:l2D]
    
    l1d,l2d = l12(7,ncoo,nbo)  
    jointpbod.rho .= modpar[l1d:l2d]
    
    return jointpbod
end
=#
##########################################

#function vecmodpar2jointstruct(jointprob::Joint2DpolyProb,modpar::Vector{<:Real})
"""
$(TYPEDSIGNATURES)

Function to reconstruct a `JointPolygBodies2D` polygons structure from i) both `Mag2DPolyMisf` and `Grav2DPolyMisf` misfit structures, ii) a list of `bodyindices` and iii) a vector of model parameters the user would like to invert for.  
"""
function vecmodpar2jointstruct(magmisf::Mag2DPolyMisf,gravmisf::Grav2DPolyMisf,
                               curbodyindices::Vector{<:Vector{<:Integer}},modpar::Vector{<:Real})
    ## separate "curbodyindices" is necessary because the bodyindices might not be in sync
    ##  with the iterations of HMC, so they must be taken from the .h5 file.

    whichparmag = magmisf.whichpar
    whichpargrav = gravmisf.whichpar
    @assert magmisf.bodyindices==gravmisf.bodyindices

    nbo = length(magmisf.bodyindices)

    if whichparmag==:all && whichpargrav==:all
        
        ncoo = length(modpar)-(7*nbo)
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

        # set density
        l1d,l2d = l12(7,ncoo,nbo)
        rho = modpar[l1d:l2d]

    elseif whichparmag==:all && whichpargrav==:vertices

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

        # set density
        rho = copy(gravmisf.rho)

    elseif whichparmag==:vertices && whichpargrav==:all
        
        ncoo = length(modpar)-nbo
        nvert = div(ncoo,2)

        # set vertices
        allvert = reshape(modpar[1:ncoo],nvert,2)

        # set magnetization
        Jind = deepcopy(magmisf.Jind)
        Jrem = deepcopy(magmisf.Jrem)

        # set density
        rho = modpar[ncoo+1:end]
        
    elseif whichparmag==:magnetization && whichpargrav==:density

        # set vertices
        @assert magmisf.allvert==gravmisf.allvert
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

        # set density
        rho = modpar[(6*nbo)+1:end]
        
    elseif whichparmag==:vertices && whichpargrav==:vertices

        ncoo = length(modpar)
        nvert = div(ncoo,2)

        # set vvertices
        allvert = reshape(modpar,nvert,2)

        # set magnetization
        Jind = deepcopy(magmisf.Jind)
        Jrem = deepcopy(magmisf.Jrem)

        # set density
        rho = copy(gravmisf.rho)

    else
        error("vecmodpar2jointstruct(): Wrong argument 'magmisf.whichpar' and/or 'gravmisf.whichpar': possible choices are: 
                          - magmisf.whichpar==:vertices && gravmisf.whichpar==:vertices;
                          - magmisf.whichpar==:vertices && gravmisf.whichpar==:all;
                          - magmisf.whichpar==:all && gravmisf.whichpar==:all;
                          - magmisf.whichpar==:all && gravmisf.whichpar==:vertices;
                          - magmisf.whichpar==:magnetization && gravmisf.whichpar==:density.
                          Aborting!")
    end
    
    mbody = JointPolygBodies2D(curbodyindices,allvert,Jind,Jrem,rho,
                               ylatext=magmisf.ylatext)

    return mbody
end

#############################################
##########################################

"""
$(TYPEDSIGNATURES)

Function to extract a list of polygon vertices from i) both `Mag2DPolyMisf` and `Grav2DPolyMisf` misfit structures and ii) a vector of model parameters the user would like to invert for.
"""
function vecmodpar2jointvertices(magmisf::Mag2DPolyMisf,gravmisf::Grav2DPolyMisf,modpar::Vector{<:Real})

    whichparmag = magmisf.whichpar
    whichpargrav = gravmisf.whichpar
    @assert magmisf.bodyindices==gravmisf.bodyindices

    nbo = length(magmisf.bodyindices)

    if whichparmag==:all && whichpargrav==:all
        
        ncoo = length(modpar)-(7*nbo)
        nvert = div(ncoo,2)

        # set vertices
        allvert = reshape(modpar[1:ncoo],nvert,2)

    elseif whichparmag==:all && whichpargrav==:vertices

        ncoo = length(modpar)-(6*nbo)
        nvert = div(ncoo,2)

        # set vertices
        allvert = reshape(modpar[1:ncoo],nvert,2)

    elseif whichparmag==:vertices && whichpargrav==:all
        
        ncoo = length(modpar)-nbo
        nvert = div(ncoo,2)

        # set vertices
        allvert = reshape(modpar[1:ncoo],nvert,2)
        
    elseif whichparmag==:magnetization && whichpargrav==:density

        # set vertices
        @assert magmisf.allvert==gravmisf.allvert
        allvert = copy(magmisf.allvert)

    elseif whichparmag==:vertices && whichpargrav==:vertices

        ncoo = length(modpar)
        nvert = div(ncoo,2)

        # set vvertices
        allvert = reshape(modpar,nvert,2)

    else
        error("vecmodpar2jointvertices(): Wrong argument 'magmisf.whichpar' and/or 'gravmisf.whichpar'. Aborting.")
    end
    

    return allvert
end

############################################################################
