
#######################################################

function splitmaggrav(magmisf::Mag2DPolyMisf,gravmisf::Grav2DPolyMisf,vecmodpar::Vector)

    nbi = length(magmisf.bodyindices)
    @assert nbi==length(gravmisf.bodyindices)

    if magmisf.whichpar==:all && gravmisf.whichpar==:all

        # mag
        vecmodmag = copy(vecmodpar[1:end-nbi])
        # grav
        vecmodgrav = vcat(vecmodpar[1:end-(7*nbi)],vecmodpar[end-nbi+1:end])
            
    elseif magmisf.whichpar==:all && gravmisf.whichpar==:vertices

        # mag
        vecmodmag = copy(vecmodpar[1:end])
        # grav
        vecmodgrav = copy(vecmodpar[1:end-(6*nbi)])

    elseif magmisf.whichpar==:vertices && gravmisf.whichpar==:all

        # mag
        vecmodmag = copy(vecmodpar[1:end-nbi])
        # grav
        vecmodgrav = copy(vecmodpar[1:end])

    elseif magmisf.whichpar==:magnetization && gravmisf.whichpar==:density
        
        # mag
        vecmodmag = copy(vecmodpar[1:end-nbi])
        # grav
        vecmodgrav = copy(vecmodpar[end-nbi+1:end])
        
    elseif magmisf.whichpar==:vertices && gravmisf.whichpar==:vertices

        # mag
        vecmodmag = copy(vecmodpar) 
        # grav
        vecmodgrav = copy(vecmodpar)
        
    end
    
    return vecmodmag,vecmodgrav
end

#########################################################################

"""
$(TYPEDSIGNATURES)

 Calculate the misfit value.
"""
function calcmisfjointmaggrav(magmisf::Mag2DPolyMisf,gravmisf::Grav2DPolyMisf,
                              vecmodpar::AbstractArray)

    vecmodmag,vecmodgrav = splitmaggrav(magmisf,gravmisf,vecmodpar)
    
    misvalmag  = magmisf(vecmodmag)
    misvalgrav = gravmisf(vecmodgrav)
    misval = misvalmag+misvalgrav

    return misval
end

#########################################################################

"""
$(TYPEDSIGNATURES)

 Pre-calculate some parameters before computing the gradient of the misfit using automatic differentiation.
"""
function precalcADstuffjointmaggrav(magmisf::Mag2DPolyMisf,gravmisf::Grav2DPolyMisf,ADkindmag::String,
                                    ADkindgrav::String,vecmodpar::AbstractArray)
    
    vecmodmag,vecmodgrav = splitmaggrav(magmisf,gravmisf,vecmodpar)

    autodiffstuffmag = MagGrav2Dpoly.precalcADstuffmag(magmisf,ADkindmag,vecmodmag)

    autodiffstuffgrav = MagGrav2Dpoly.precalcADstuffgrav(gravmisf,ADkindgrav,vecmodgrav)
    
    return autodiffstuffmag,autodiffstuffgrav
end   

#######################################################
"""
$(TYPEDSIGNATURES)

Function to compute the gradient of the misfit with respect to the model parameters as required for HMC inversions. 
The gradient is computed by means of automatic differentiation using one of three different methods. The user must indicate the method by `ADkind` string, choosing among `FWDdiff`, `REVdiffTAPE` or `REVdiffTAPEcomp`. 
For an explanation about the automatic differentiation method, the reader is invited to look at the documentation 
relative to the Julia packages `ForwardDiff` and `ReverseDiff`.  
"""  
function calc∇misfjointmaggrav(magmisf::Mag2DPolyMisf,gravmisf::Grav2DPolyMisf,ADkindmag::String,
                               ADkindgrav::String,vecmodpar::AbstractArray,
                               autodiffstuffmag,autodiffstuffgrav)

    vecmodmag,vecmodgrav = splitmaggrav(magmisf,gravmisf,vecmodpar)
    
    gradmag = calc∇misfmag(magmisf,vecmodmag,ADkindmag,autodiffstuffmag)

    gradgrav = calc∇misfgrav(gravmisf,vecmodgrav,ADkindgrav,autodiffstuffgrav)

    # merge mag & grav gradients
    nbi = length(magmisf.bodyindices)
    @assert nbi==length(gravmisf.bodyindices)


    if magmisf.whichpar==:all && gravmisf.whichpar==:all
        
        gradver = gradmag[1:end-(6*nbi)]+gradgrav[1:end-nbi]
        gradmgn = gradmag[end-(6*nbi)+1:end]
        graddens= gradgrav[end-nbi+1:end]
        grad = vcat(gradver,gradmgn,graddens)
        
    elseif magmisf.whichpar==:all && gravmisf.whichpar==:vertices

        gradver = gradmag[1:end-(6*nbi)]+gradgrav[1:end]
        gradmgn = gradmag[end-(6*nbi)+1:end]
        grad = vcat(gradver,gradmgn)
        
    elseif magmisf.whichpar==:vertices && gravmisf.whichpar==:all

        gradver = gradmag[1:end]+gradgrav[1:end-nbi]
        graddens = gradgrav[end-nbi+1:end]
        grad = vcat(gradver,graddens)
        
    elseif magmisf.whichpar==:magnetization && gravmisf.whichpar==:density
        
        grad = vcat(gradmag,gradgrav)
        
    elseif magmisf.whichpar==:vertices && gravmisf.whichpar==:vertices
        
        grad = gradmag+gradgrav
        
    end

    
    nangrad = isnan.(grad)
    if any(nangrad)
        error("calc∇misfjointmaggrav(): isnan.(grad) = $nangrad")
    end
    return grad
end

#######################################################
