
"""
HMCGrav2Dpoly

A convenience module to facilitate the use of `Grav2Dpoly` within the framework of Hamiltonian Monte Carlo inversion by employing the package `HMCtomo`. 

# Exports

$(EXPORTS)
"""
module HMCGrav2Dpoly

using MagGrav2Dpoly
using GeoPolygons
using ForwardDiff
using ReverseDiff
using DocStringExtensions
using HDF5
using LinearAlgebra: LowerTriangular,Diagonal
#using Random

#using Diffractor

export Grav2DpolyProb,gravprob2vec

#################################################################

## create the problem type
"""
$(TYPEDSIGNATURES)

Structure to define a 'problem' for HMC. 
The users must indicate the method to compute the misfit gradient by setting the variable `ADkind`, choosing among `FWDdiff`, `REVdiffTAPE` or `REVdiffTAPEcomp`. 
For an explanation of the automatic differentiation method, the reader is invited to look at the documentation 
relative to the Julia packages `ForwardDiff` and `ReverseDiff`. 
"""
struct Grav2DpolyProb
    gravmisf::Grav2DPolyMisf
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
    
    #function Grav2DpolyProb(gravmisf,topography,M,ADkind)
    function Grav2DpolyProb(; pbodystart::GravPolygBodies2D,topography::TopoEdges,
                            grav_xzobs::Array{<:Real,2},grav_obsdata::Vector{<:Real},
                            grav_invCd::AbstractMatrix{<:Real},grav_whichpar::Symbol,grav_ADkind::String,
                            trytofixpolygons::Bool=false )

        # Create  Grav2DPolyMisf
        if grav_whichpar==:all
            grav_allvert=nothing
            rho=nothing
        elseif grav_whichpar==:vertices
            grav_allvert=nothing
            rho=copy(pbodystart.rho)
        elseif grav_whichpar==:density
            grav_allvert=deepcopy(pbodystart.geom.allvert)
            rho=nothing
        end
        
        if pbodystart.ylatext==nothing
            ylatext = nothing
        else
            ylatext = copy(pbodystart.ylatext)
        end
        
        if ylatext==nothing
            printstyled(" Gravity problem type: 2D\n",bold=false,color=:yellow)
        elseif -ylatext[1]==ylatext[2]
            printstyled(" Gravity problem type: 2.5D\n",bold=false,color=:yellow)
        else
            printstyled(" Gravity problem type: 2.75D\n",bold=false,color=:yellow)
        end

        # Instantiate the misfit type
        gravmisf = Grav2DPolyMisf(pbodystart.geom.bodyindices,grav_xzobs,grav_obsdata,grav_invCd,
                                  grav_whichpar,allvert=grav_allvert,rho=rho,ylatext=ylatext)

        ##----------------------------------
        firsttime_grad = Ref(true) 
        autodiffstuff  = Ref(nothing)  ## This fix should be improved!!!
        dofwdchecks = true
        firstcheck = Ref(true)
        orgbodyindices = deepcopy(gravmisf.bodyindices)
        probkind = :polygonalbodies
        #--------------------------------
        # Print whichpar for grav misf
        if gravmisf.whichpar==:all
            println(" ---------------------------------------------------")
            printstyled("\n Type Parameters to Invert: Vertices + Density\n",bold=true,color=:light_red)
            println("\n ---------------------------------------------------")
        elseif gravmisf.whichpar==:vertices
            println(" ---------------------------------------------------")
            printstyled("\n Type Parameters to Invert: Vertices\n",bold=true,color=:light_red)
            println("\n ---------------------------------------------------")
        elseif gravmisf.whichpar==:density
            println(" ---------------------------------------------------")
            printstyled("\n Type Parameters to Invert: Density\n",bold=true,color=:light_red)
            println("\n ---------------------------------------------------")
        end
        
        return new(gravmisf,topography,grav_ADkind,firsttime_grad,autodiffstuff,
                   dofwdchecks,firstcheck,orgbodyindices,probkind,trytofixpolygons)
    end
end

## use  x.T * C^-1 * x  = ||L^-1 * x ||^2 ?
#################################################################

## make the type callable
function (grav2dprob::Grav2DpolyProb)(vecmodpar::Vector{Float64},kind::Symbol)

 
    if kind==:nlogpdf
        #############################################
        ## compute the logdensity value for vecvel ##
        #############################################
        #println("logpdf")
        misval = grav2dprob.gravmisf(vecmodpar)
        return misval
        

    elseif kind==:gradnlogpdf
        #################################################
        ## compute the gradient of the misfit function ##
        #################################################

        ## create the REVdiff tape, compile and store it using a pointer
        #usingrevdiff = Bool(grav2dprob.ADkind=="REVdiffTAPE" || grav2dprob.ADkind=="REVdiffTAPEcomp")

        if grav2dprob.firsttime_grad[] #&& (usingrevdiff==true)

            grav2dprob.autodiffstuff[] = Grav2Dpoly.precalcADstuffgrav(gravmisf,ADkind,vecmodpar)

            # if grav2dprob.ADkind=="REVdiffTAPE"
            #     grav2dprob.autodiffstuff[] = ReverseDiff.GradientTape(grav2dprob.gravmisf,vecmodpar)

            # elseif grav2dprob.ADkind=="REVdiffTAPEcomp"
            #     # compiled tape
            #     tape = ReverseDiff.GradientTape(grav2dprob.gravmisf,vecmodpar)
            #     grav2dprob.autodiffstuff[] = ReverseDiff.compile(tape)

            # elseif grav2dprob.ADkind=="FWDdiff"
            #     if length(vecmodpar)<20
            #         chuncksize = length(vecmodpar)
            #     else
            #         chuncksize = 20 # 20
            #     end
            #     cfggrd = ForwardDiff.GradientConfig(grav2dprob.gravmisf,vecmodpar,ForwardDiff.Chunk{chuncksize}())
            #     grav2dprob.autodiffstuff[] = cfggrd

            # else
            #     grav2dprob.autodiffstuff[] = nothing
                
            # end
            
            ## set first time to false!
            grav2dprob.firsttime_grad[]=false
        end

        # compute gradient
        vecgrad = Grav2Dpoly.calc∇misfgrav(grav2dprob.gravmisf,vecmodpar,
                                           grav2dprob.ADkind,grav2dprob.autodiffstuff[])
        
        # return flattened gradient
        return vecgrad
        
    else
        error("grav2dprob::Grav2DpolyProb(): Wrong argument 'kind': $kind...")
    end
end

#################################################################

"""
 Perform tests on polygons for Grav2D HMC inversion.
"""
function (gravprob::Grav2DpolyProb)(mcur::Vector{Float64},mnew::Vector{Float64},pnew::Vector{Float64},
                                    lowcon::Union{Vector{Float64},Nothing},upcon::Union{Vector{Float64},Nothing},
                                    LcholmassM::Union{LowerTriangular{Float64,Matrix{Float64}},Diagonal{Float64}},
                                    dowhat::String)

    if dowhat=="performchecks"

        # perform checks
        allgood = hmcpolycheckgrav!(gravprob,mcur,mnew,pnew,lowcon,upcon,LcholmassM)

    end
    
    return allgood
end

#################################################################

function (gravprob::Grav2DpolyProb)(dowhat::String;
                                    newbodyindices::Union{Nothing,Vector{<:Vector{<:Integer}}}=nothing)

    if dowhat=="newinds"
        if newbodyindices==nothing
            gravprob.orgbodyindices .= gravprob.gravmisf.bodyindices
        else
            gravprob.orgbodyindices .= newbodyindices
        end
        
    elseif dowhat=="backinds"
        gravprob.gravmisf.bodyindices .= gravprob.orgbodyindices

    elseif dowhat=="getbodyindices"
        # COPY the bodyindices otherwise they'll get modified!
        return copy(gravprob.gravmisf.bodyindices) 
    end
    
    return 
end

######################################################

function (gravprob::Grav2DpolyProb)(fid_h5::HDF5.File,dowhat::String)

    if dowhat=="readh5bodyindices"

        nbo = length(gravprob.orgbodyindices)
        for i=1:nbo
            dset = "curbodyindices_$i"
            gravprob.orgbodyindices[i] = fid_h5[dset][:,end]
            gravprob.gravmisf.bodyindices[i] = fid_h5[dset][:,end]
        end

    end
    return 
end

##############################################################

function (gravprob::Grav2DpolyProb)(fid_h5::HDF5.File,it::Integer,sit::Integer,dowhat::String)

    if dowhat=="savebodyindices"

        nbo = length(gravprob.orgbodyindices)
        if it==1
            for i=1:nbo
                nind = length(gravprob.orgbodyindices[i])
                bi_h5 = HDF5.create_dataset(fid_h5,"curbodyindices_$i", Int64, ((nind,1), (nind,-1)),chunk=(nind,1))
                HDF5.h5d_close(bi_h5)
            end
        end
    
        # models
        for i=1:nbo
            
            nind = length(gravprob.orgbodyindices[i])
            bi_h5 = HDF5.open_dataset(fid_h5,"curbodyindices_$i")
            HDF5.set_extent_dims(bi_h5,(nind,sit))
            bi_h5[:,sit] = gravprob.orgbodyindices[i][:]
            HDF5.h5d_close(bi_h5)

        end

    end
    return 
end

##############################################################

function gravprob2vec(gravprob::Grav2DpolyProb,gravpbod::GravPolygBodies2D)
    gravstruct2vec(gravprob.gravmisf.whichpar,gravpbod)
end

#########################################################################

"""
$(TYPEDSIGNATURES)

Function carrying out all the checks in the case of model parameterization 
as polygonal bodies in the framework of HMC inversion.     
"""
function hmcpolycheckgrav!(gravprob::Grav2DpolyProb,mcur::Vector{Float64},mnew::Vector{Float64},pnew::Vector{Float64},
                           lowcon::Union{Vector{Float64},Nothing},upcon::Union{Vector{Float64},Nothing},
                           LcholmassM::Union{LowerTriangular{Float64,Matrix{Float64}},Diagonal{Float64}};
                           maxclk::Integer=3)

    ##===========================================================
    ## CHECKS on densities, vertices etc. in case of polygons
    ##===========================================================
    # No checks if whichpar = density and initial checks already performed
    if !gravprob.firstcheck[] && gravprob.gravmisf.whichpar==:density
        return true
    end

    #-------------------------
    ## create some structures
    qbody = Grav2Dpoly.vecmodpar2gravstruct(gravprob.gravmisf,mnew)
    
    ##-------------------------------------------------
    # Initial check intersections of segments and with topography
    #--> Aborting in ase of intersections
    if gravprob.firstcheck[]
        intersect = GeoPolygons.checkall(qbody.geom.bo,gravprob.topography)
        if !all(intersect.==:none)
            error("Possible either interesections between polygons or crossings with topography in your starting model. Re-check your parameterization! Aborting!")
        end
        gravprob.firstcheck[]=false
        return nothing
    end
    
    ## fix problems with topography and intersections
    clk = 0
        
    ##-------------------------------------------------
    ## check ordering of vertices (must be anticlockwise)
    for k=1:qbody.geom.nbo
        clockw = checkanticlockwiseorder(qbody.geom.bo[k])
        
        if !clockw
            gravprob.gravmisf.bodyindices[k] = reverse(gravprob.gravmisf.bodyindices[k],dims = 1)
        end
    end
    
    qbody = GravPolygBodies2D(gravprob.gravmisf.bodyindices,qbody.geom.allvert,qbody.rho,
                              ylatext=gravprob.gravmisf.ylatext)
    
    #----------------------------------------------------
    # check intersections of segments and with topography
    intersect = GeoPolygons.checkall(qbody.geom.bo,gravprob.topography)
    
    if gravprob.trytofixpolygons

        ### attempt to fix polygons
        while !all(intersect.==:none)
            
            # performing all fixing in a random order
            GeoPolygons.fixall!(gravprob.topography,qbody.geom.bo,intersect,gravprob.gravmisf.bodyindices,lowcon,upcon)
            
            ##-------------------------------------------------
            ## check ordering of vertices (must be anticlockwise)
            for k=1:qbody.geom.nbo
                clockw = checkanticlockwiseorder(qbody.geom.bo[k])
                
                if !clockw
                    gravprob.gravmisf.bodyindices[k] = reverse(gravprob.gravmisf.bodyindices[k],dims = 1)
                end
            end
            
            qbody = GravPolygBodies2D(gravprob.gravmisf.bodyindices,qbody.geom.allvert,qbody.rho,
                                      ylatext=gravprob.gravmisf.ylatext)        
            # do intersection checks
            intersect = GeoPolygons.checkall(qbody.geom.bo,gravprob.topography)
            
            clk += 1
            
            if clk > maxclk
                # check succeded? nope...
                gravprob("backinds")
                return false
            end    
        end

    else

        ## do NOT attempt to fix polygons        
        if !all(intersect.==:none)
            # check succeded? nope...
            gravprob("backinds")
            return false
        end

    end
    
    ##-------------------------------------------------
    ## unroll back to vectors
    mnew .= Grav2Dpoly.gravstruct2vec(gravprob.gravmisf.whichpar,qbody)
    
    ## Adjust p to the new altered trajectory
    ## pnew = M * grad(K) [ from gradk(K) = M^-1 * pnew ]
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
