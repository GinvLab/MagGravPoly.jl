
"""
HMCMag2Dpoly

A convenience module to facilitate the use of the magnetic forward code in `MG2D` within the framework of Hamiltonian Monte Carlo inversion by employing the package `MCsamplers`. 

# Exports

$(EXPORTS)
"""
module HMCMag2Dpoly

using ..MG2D
using ..GeoPoly
using ForwardDiff
using ReverseDiff
using DocStringExtensions
using HDF5
using LinearAlgebra: LowerTriangular,Diagonal
#using Random  

#using Diffractor

export Mag2DpolyProb,magprob2vec


#################################################################

## create the problem type
"""
$(TYPEDSIGNATURES)

Julia structure to define a magnetic problem for HMC.
The users must indicate the method to compute misfit gradient by `ADkind` string, choosing among `FWDdiff`, `REVdiffTAPE` or `REVdiffTAPEcomp`. 
For an explanation of the automatic differentiation method, the reader is invited to look at the documentation 
relative to the Julia packages `ForwardDiff` and `ReverseDiff`. 
"""
struct Mag2DpolyProb
    magmisf::Mag2DPolyMisf
    topography::TopoEdges
    ADkind::String
    firsttime_grad::Base.RefValue{Bool} 
    autodiffstuff::Base.Ref{Any}  ## This fix should be improved!!!
    dofwdchecks::Bool  # perform checks during HMC leapfrog steps
    firstcheck::Base.RefValue{Bool} 
    orgbodyindices::Vector{Vector{<:Integer}}
    #M::AbstractMatrix{Float64}
    probkind::Symbol
    trytofixpolygons::Bool
    
    ##-------------------------------------------------------------
    #function Mag2DpolyProb(magmisf,topography,M,ADkind)
    function Mag2DpolyProb(; pbodystart::MagPolygBodies2D,topography::TopoEdges,
                           northxax::Real,mag_xzobs::Array{<:Real,2},mag_obsdata::Vector{<:Real},
                           mag_invCd::AbstractMatrix{<:Real},mag_whichpar::Symbol,mag_ADkind::String, 
                           trytofixpolygons::Bool=false )

        # First create Mag2DPolyMisf 
        if mag_whichpar==:all
            mag_allvert=nothing
            Jind = nothing
            Jrem = nothing
        elseif mag_whichpar==:vertices
            mag_allvert=nothing
            Jind = deepcopy(pbodystart.Jind)
            Jrem = deepcopy(pbodystart.Jrem)
        elseif mag_whichpar==:magnetization
            mag_allvert = copy(pbodystart.geom.allvert)
            Jind = nothing
            Jrem = nothing
        end
        if pbodystart.ylatext==nothing
            ylatext = nothing
        else
            ylatext = copy(pbodystart.ylatext)
        end

        println(" ---------------------------------------------------")
        if ylatext==nothing
            printstyled(" Magnetic problem type: 2D\n",bold=true,color=:light_yellow)
        elseif -ylatext[1]==ylatext[2]
            printstyled(" Magnetic problem type: 2.5D\n",bold=true,color=:light_yellow)
        else
            printstyled(" Magnetic problem type: 2.75D\n",bold=true,color=:light_yellow)
        end

        # Instantiate the misfit type
        magmisf = Mag2DPolyMisf(pbodystart.geom.bodyindices,northxax,mag_xzobs,mag_obsdata,mag_invCd,mag_whichpar,
                                allvert=mag_allvert,Jind=Jind,Jrem=Jrem,ylatext=ylatext)

        ##----------------------------------------------
        firsttime_grad = Ref(true) 
        autodiffstuff  = Ref(nothing)  ## This fix should be improved!!!
        dofwdchecks = true
        firstcheck = Ref(true)
        orgbodyindices = deepcopy(magmisf.bodyindices)
        probkind = :polygonalbodies
        #--------------------------------
        # Print whichpar for mag misf
        if magmisf.whichpar==:all
            println(" ---------------------------------------------------")
            printstyled(" Parameters to invert for: Vertices + Magnetization\n",bold=true,color=:light_red)
        elseif magmisf.whichpar==:vertices
            println(" ---------------------------------------------------")
            printstyled(" Parameters to invert for: Vertices\n",bold=true,color=:light_red)
        elseif magmisf.whichpar==:magnetization
            println(" ---------------------------------------------------")
            printstyled(" Parameters to invert for: Magnetization\n",bold=true,color=:light_red)
        end
        
        return new(magmisf,topography,mag_ADkind,firsttime_grad,autodiffstuff,
                   dofwdchecks,firstcheck,orgbodyindices,probkind,trytofixpolygons)
    end
end

## use  x.T * C^-1 * x  = ||L^-1 * x ||^2 ?
###########################################################################

## make the type callable
function (mag2dprob::Mag2DpolyProb)(vecmodpar::Vector{Float64},kind::Symbol)

 
    if kind==:nlogpdf
        #############################################
        ## compute the logdensity value for vecvel ##
        #############################################
        #println("logpdf")
        misval = mag2dprob.magmisf(vecmodpar)
        return misval
        

    elseif kind==:gradnlogpdf
        #################################################
        ## compute the gradient of the misfit function ##
        #################################################

        ## create the REVdiff tape, compile and store it using a pointer
        #usingrevdiff = Bool(mag2dprob.ADkind=="REVdiffTAPE" || mag2dprob.ADkind=="REVdiffTAPEcomp")

        if mag2dprob.firsttime_grad[] #&& (usingrevdiff==true)

            mag2dprob.autodiffstuff[] = MG2D.precalcADstuffmag(mag2dprob.magmisf,mag2dprob.ADkind,vecmodpar)

            # if mag2dprob.ADkind=="REVdiffTAPE"
            #     mag2dprob.autodiffstuff[] = ReverseDiff.GradientTape(mag2dprob.magmisf,vecmodpar)

            # elseif mag2dprob.ADkind=="REVdiffTAPEcomp"
            #     # compiled tape
            #     tape = ReverseDiff.GradientTape(mag2dprob.magmisf,vecmodpar)
            #     mag2dprob.autodiffstuff[] = ReverseDiff.compile(tape)

            # elseif mag2dprob.ADkind=="FWDdiff"
            #     if length(vecmodpar)<20
            #         chuncksize = length(vecmodpar)
            #     else
            #         chuncksize = 20 # 20
            #     end
            #     cfggrd = ForwardDiff.GradientConfig(mag2dprob.magmisf,vecmodpar,ForwardDiff.Chunk{chuncksize}())
            #     mag2dprob.autodiffstuff[] = cfggrd

            # else
            #     mag2dprob.autodiffstuff[] = nothing
                
            # end
            
            ## set first time to false!
            mag2dprob.firsttime_grad[]=false
        end

        # compute gradient
        vecgrad = MG2D.calc∇misfmag(mag2dprob.magmisf,vecmodpar,
                        mag2dprob.ADkind,mag2dprob.autodiffstuff[])
        
        # return flattened gradient
        return vecgrad
        
    else
        error("mag2dprob::Mag2DpolyProb(): Wrong argument 'kind': $kind...")
    end
end

#################################################################

"""
 Perform tests on polygons for Mag2D HMC inversion.
"""
function (magprob::Mag2DpolyProb)(mcur::Vector{Float64},mnew::Vector{Float64},pnew::Vector{Float64},
                                  lowcon::Union{Vector{Float64},Nothing},upcon::Union{Vector{Float64},Nothing},
                                  LcholmassM::Union{LowerTriangular{Float64,Matrix{Float64}},Diagonal{Float64}},
                                  dowhat::String)

    if dowhat=="performchecks"

        # perform checks
        allgood = hmcpolycheckmag!(magprob,mcur,mnew,pnew,lowcon,upcon,LcholmassM)

    end
    
    return allgood
end

##############################################################

function (magprob::Mag2DpolyProb)(dowhat::String ;
                                  newbodyindices::Union{Nothing,Vector{<:Vector{<:Integer}}}=nothing)
    
    if dowhat=="newinds"
        if newbodyindices==nothing
            magprob.orgbodyindices .= magprob.magmisf.bodyindices
        else
            magprob.orgbodyindices .= newbodyindices
        end

    elseif dowhat=="backinds"
        magprob.magmisf.bodyindices .= magprob.orgbodyindices

    elseif dowhat=="getbodyindices"
        # COPY the bodyindices otherwise they'll get modified!
        return copy(magprob.magmisf.bodyindices) 
    end

    return
end

######################################################

function (magprob::Mag2DpolyProb)(fid_h5::HDF5.File,dowhat::String)

    if dowhat=="readh5bodyindices"

        nbo = length(magprob.orgbodyindices)
        for i=1:nbo
            dset = "curbodyindices_$i"
            magprob.orgbodyindices[i] = fid_h5[dset][:,end]
            magprob.magmisf.bodyindices[i] = fid_h5[dset][:,end]
        end

    end
    return 
end

##############################################################

function (magprob::Mag2DpolyProb)(fid_h5::HDF5.File,it::Integer,sit::Integer,dowhat::String)

    if dowhat=="savebodyindices"

        nbo = length(magprob.orgbodyindices)
        if it==1
            for i=1:nbo
                nind = length(magprob.orgbodyindices[i])
                bi_h5 = HDF5.create_dataset(fid_h5,"curbodyindices_$i", Int64, ((nind,1), (nind,-1)),chunk=(nind,1))
                HDF5.h5d_close(bi_h5)
            end
        end
    
        # models
        for i=1:nbo
            
            nind = length(magprob.orgbodyindices[i])
            bi_h5 = HDF5.open_dataset(fid_h5,"curbodyindices_$i")
            HDF5.set_extent_dims(bi_h5,(nind,sit))
            bi_h5[:,sit] = magprob.orgbodyindices[i][:]
            HDF5.h5d_close(bi_h5)
            
        end

    end
    return 
end

##############################################################

function magprob2vec(magprob::Mag2DpolyProb,magpbod::MagPolygBodies2D)
    magstruct2vec(magprob.magmisf.whichpar,magpbod)
end

#########################################################################

"""
$(TYPEDSIGNATURES)

Function carrying out all the checks in the case of model parameterization 
as polygonal bodies in the framework of HMC inversion.     
"""
function hmcpolycheckmag!(magprob::Mag2DpolyProb,mcur::Vector{Float64},mnew::Vector{Float64},pnew::Vector{Float64},
                          lowcon::Union{Vector{Float64},Nothing},upcon::Union{Vector{Float64},Nothing},
                          LcholmassM::Union{LowerTriangular{Float64,Matrix{Float64}},Diagonal{Float64}};
                          maxclk::Integer=3)

    ##===========================================================
    ## CHECKS on magnetization, vertices etc. in case of polygons
    ##===========================================================
    # No checks if whichpar = magnetization and initial checks already performed
    if !magprob.firstcheck[] && magprob.magmisf.whichpar==:magnetization
        return true
    end
   
    #-------------------------
    ## create some structures
    qbody = MG2D.vecmodpar2magstruct(magprob.magmisf,magprob.magmisf.bodyindices,mnew)

    ##-------------------------------------------------
    # Initial check intersections of segments and with topography
    #--> Aborting in case of intersections
    if magprob.firstcheck[]
        intersect = GeoPoly.checkall(qbody.geom.bo,magprob.topography)
        if !all(intersect.==:none)
            error("Possible either interesections between polygons or crossings with topography in your starting model. Re-check your parameterization! Aborting!")
        end
        magprob.firstcheck[]=false
        return nothing
    end
        
    ## fix problems with topography and intersections
    clk = 0
    
    ##-------------------------------------------------
    ## check ordering of vertices (must be anticlockwise)
    for k=1:qbody.geom.nbo
        clockw = checkanticlockwiseorder(qbody.geom.bo[k])
        
        if !clockw  
            magprob.magmisf.bodyindices[k] = reverse(magprob.magmisf.bodyindices[k],dims = 1)
        end
    end
    
    qbody = MagPolygBodies2D(magprob.magmisf.bodyindices,qbody.geom.allvert,qbody.Jind,qbody.Jrem,
                             ylatext=magprob.magmisf.ylatext)
    
    #----------------------------------------------------
    # check intersections of segments and with topography
    intersect = GeoPoly.checkall(qbody.geom.bo,magprob.topography)
    
    if magprob.trytofixpolygons

        ### attempt to fix polygons
        while !all(intersect.==:none)
            
            # performing all fixing in a random order
            GeoPoly.fixall!(magprob.topography,qbody.geom.bo,intersect,magprob.magmisf.bodyindices,lowcon,upcon)
            
            ##-------------------------------------------------
            ## check ordering of vertices (must be anticlockwise)
            for k=1:qbody.geom.nbo
                clockw = checkanticlockwiseorder(qbody.geom.bo[k])
                
                if !clockw
                    magprob.magmisf.bodyindices[k] = reverse(magprob.magmisf.bodyindices[k],dims = 1)
                end                
            end
            
            qbody = MagPolygBodies2D(magprob.magmisf.bodyindices,qbody.geom.allvert,qbody.Jind,qbody.Jrem,
                                     ylatext=magprob.magmisf.ylatext)
            
            # do intersection checks
            intersect = GeoPoly.checkall(qbody.geom.bo,magprob.topography)
            
            clk += 1
            
            if clk > maxclk
                # check succeded? nope...
                magprob("backinds")
                return false
            end    
        end

    else

        ## do NOT attempt to fix polygons        
        if !all(intersect.==:none)
            # check succeded? nope...
            magprob("backinds")
            return false
        end

    end
        
    ##-------------------------------------------------
    ## unroll back to vectors
    mnew .= MG2D.magstruct2vec(magprob.magmisf.whichpar,qbody)

    ## Adjust p to the new altered trajectory
    ## pnew = M * grad(K) [ from gradk(K) = M^-1 * p ]
    ##  where grad(K) = mnew - mcur
    if clk > 0
        grak = mnew .- mcur
        ## parenthesize such that there are two matrix-vector products
        pnew = LcholmassM * ( transpose(LcholmassM) * grak ) 
    end

    # check succeded? yes!
    return true
end

########################################################################################


end # module
#################################################################
