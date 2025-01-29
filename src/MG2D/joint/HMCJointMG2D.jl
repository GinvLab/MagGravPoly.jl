
"""
HMCJointMG2D

A convenience module to facilitate the use of both magnetic and gravity forward codes in `MG2D` within the framework of Hamiltonian Monte Carlo inversion by employing the package `InverseAlgos.Samplers.jl`. 

# Exports

$(EXPORTS)
"""
module HMCJointMG2D 

using ..MG2D
using ..GeoPoly
using ForwardDiff
using ReverseDiff
using DocStringExtensions
using HDF5
using LinearAlgebra: LowerTriangular,Diagonal

export Joint2DpolyProb,jointprob2vec

#################################################################

## create the problem type
"""
$(TYPEDSIGNATURES)

Julia structure to define a joint magnetic and gravity problem for HMC. 
The users must indicate the method to compute misfit gradient for both the magnetic and gravity problem by `ADkindmag` and `ADkindgrav` strings, respectively, choosing among `FWDdiff`, `REVdiffTAPE` or `REVdiffTAPEcomp`. 
For an explanation of the automatic differentiation method, the reader is invited to look at the documentation 
relative to the Julia packages `ForwardDiff` and `ReverseDiff`.
"""
struct Joint2DpolyProb
    magmisf::Mag2DPolyMisf
    gravmisf::Grav2DPolyMisf
    topography::TopoEdges
    ADkindmag::String
    ADkindgrav::String
    firsttime_grad::Base.RefValue{Bool} 
    autodiffstuffmag::Base.Ref{Any}  ## This fix should be improved!!!
    autodiffstuffgrav::Base.Ref{Any}
    dofwdchecks::Bool   # perform checks during HMC leapfrog steps
    firstcheck::Base.RefValue{Bool}
    orgbodyindices::Vector{Vector{<:Integer}}
    #M::AbstractMatrix{<:Real}
    probkind::Symbol
    trytofixpolygons::Bool

    ##-------------------------------------------------------------
    function Joint2DpolyProb(; jpbodystart::JointPolygBodies2D,topography::TopoEdges,
                             northxax::Real,mag_xzobs::Array{<:Real,2},mag_obsdata::Vector{<:Real},
                             mag_invCd::AbstractMatrix{<:Real},mag_whichpar::Symbol,mag_ADkind::String,
                             grav_xzobs::Array{<:Real,2},grav_obsdata::Vector{<:Real},
                             grav_invCd::AbstractMatrix{<:Real},grav_whichpar::Symbol,grav_ADkind::String,
                             trytofixpolygons::Bool=false )



        # First create Mag2DPolyMisf 
        if mag_whichpar==:all
            mag_allvert=nothing
            Jind=nothing
            Jrem=nothing
        elseif mag_whichpar==:vertices
            mag_allvert=nothing
            Jind=jpbodystart.Jind
            Jrem=jpbodystart.Jrem
        elseif mag_whichpar==:magnetization
            mag_allvert=jpbodystart.geom.allvert
            Jind=nothing
            Jrem=nothing
        end
        if jpbodystart.ylatext==nothing
            ylatext = nothing
        else
            ylatext = copy(jpbodystart.ylatext)
        end

        println(" ---------------------------------------------------") 
        if ylatext==nothing
            printstyled(" Joint magnetic & gravity problem type: 2D\n",bold=true,color=:light_yellow)
        elseif -ylatext[1]==ylatext[2]
            printstyled(" Joint magnetic & gravity problem type: 2.5D\n",bold=true,color=:light_yellow)
        else
            printstyled(" Joint magnetic & gravity problem type: 2.75D\n",bold=true,color=:light_yellow)
        end
        
        # Instantiate the misfit type
        magmisf = Mag2DPolyMisf(jpbodystart.geom.bodyindices,northxax,mag_xzobs,mag_obsdata,mag_invCd,mag_whichpar,
                                allvert=mag_allvert,Jind=Jind,Jrem=Jrem,
                                ylatext=jpbodystart.ylatext)

        # Create  Grav2DPolyMisf
        if grav_whichpar==:all
            grav_allvert=nothing
            rho=nothing
        elseif grav_whichpar==:vertices
            grav_allvert=nothing
            rho=jpbodystart.rho 
        elseif grav_whichpar==:density
            grav_allvert=jpbodystart.geom.allvert
            rho=nothing
        end
        
        # Instantiate the misfit type
        gravmisf = Grav2DPolyMisf(jpbodystart.geom.bodyindices,grav_xzobs,grav_obsdata,grav_invCd,
                                  grav_whichpar,allvert=grav_allvert,rho=rho,
                                  ylatext=jpbodystart.ylatext)

        ##
        firsttime_grad = Ref(true)
        autodiffstuffmag = Ref(nothing) 
        autodiffstuffgrav = Ref(nothing)
        dofwdchecks = true
        firstcheck = Ref(true)
        # bodyindices stuff
        #nbi = length(magmisf.bodyindices)
        #@assert nbi==length(gravmisf.bodyindices)
        #for i=1:nbi
        #    @assert magmisf.bodyindices[i] == gravmisf.bodyindices[i]
        #end
        orgbodyindices = deepcopy(jpbodystart.geom.bodyindices)
        # check for half dim (same for both mag & grav)
        @assert magmisf.ylatext == gravmisf.ylatext
        ##==============================================
        probkind = :polygonalbodies
        #----------------------------------------------
        # checks on whichpar for both grav and mag misf
        if magmisf.whichpar==:all && gravmisf.whichpar==:all
            println(" ---------------------------------------------------")
            printstyled(" Parameters to invert for: Vertices + Magnetization + Density\n",bold=true,color=:light_red)
        elseif magmisf.whichpar==:vertices && gravmisf.whichpar==:all
            println(" ---------------------------------------------------") 
            printstyled(" Parameters to invert for: Vertices + Density\n",bold=true,color=:light_red)
        elseif magmisf.whichpar==:all && gravmisf.whichpar==:vertices
            println(" ---------------------------------------------------")
            printstyled(" Parameters to invert for: Vertices + Magnetization\n",bold=true,color=:light_red)
        elseif magmisf.whichpar==:magnetization && gravmisf.whichpar==:density
            @assert magmisf.allvert==gravmisf.allvert
            println(" ---------------------------------------------------")
            printstyled(" Parameters to invert for: Magnetization + Density\n",bold=true,color=:light_red)
        elseif magmisf.whichpar==:vertices && gravmisf.whichpar==:vertices
            println(" ---------------------------------------------------")
            printstyled(" Parameters to invert for: Vertices\n",bold=true,color=:light_red)
        else
            error("Joint2DpolyProb(): Wrong argument 'magmisf.whichpar' and/or 'gravmisf.whichpar': possible choices are: 
                - magmisf.whichpar==:vertices && gravmisf.whichpar==:vertices;
                - magmisf.whichpar==:vertices && gravmisf.whichpar==:all;
                - magmisf.whichpar==:all && gravmisf.whichpar==:all;
                - magmisf.whichpar==:all && gravmisf.whichpar==:vertices;
                - magmisf.whichpar==:magnetization && gravmisf.whichpar==:density.
                Aborting!")
        end
    
        return new(magmisf,gravmisf,topography,mag_ADkind,grav_ADkind,firsttime_grad,
                   autodiffstuffmag,autodiffstuffgrav,dofwdchecks,firstcheck,orgbodyindices,
                   probkind,trytofixpolygons)
    end
end

###########################################################

## make the type callable
function (joint2dprob::Joint2DpolyProb)(vecmodpar::Vector{Float64},kind::Symbol)

    # split mag and grav vectors of model parameters
    vecmodmag,vecmodgrav = splitmaggrav(joint2dprob,vecmodpar)
    
    if kind==:nlogpdf
        #############################################
        ## compute the logdensity value for vecvel ##
        #############################################
        #println("logpdf")
        misvalmag  = joint2dprob.magmisf(vecmodmag)
        misvalgrav = joint2dprob.gravmisf(vecmodgrav)
        misval = misvalmag+misvalgrav
        return misval
        

    elseif kind==:gradnlogpdf
        #################################################
        ## compute the gradient of the misfit function ##
        #################################################

        ## create the REVdiff tape, compile and store it using a pointer
        #usingrevdiff = Bool(mag2dprob.ADkind=="REVdiffTAPE" || mag2dprob.ADkind=="REVdiffTAPEcomp")

        if joint2dprob.firsttime_grad[] #&& (usingrevdiff==true)

            ############### Mag ################

            joint2dprob.autodiffstuffmag[] = MG2D.precalcADstuffmag(joint2dprob.magmisf,joint2dprob.ADkindmag,
                                                                         vecmodmag)

            # if joint2dprob.ADkindmag=="REVdiffTAPE"
            #     joint2dprob.autodiffstuffmag[] = ReverseDiff.GradientTape(joint2dprob.magmisf,
            #                                                               vecmodmag)
                
            # elseif joint2dprob.ADkindmag=="REVdiffTAPEcomp"
            #     # compiled tape
            #     tapemag = ReverseDiff.GradientTape(joint2dprob.magmisf,vecmodmag)
            #     joint2dprob.autodiffstuffmag[] = ReverseDiff.compile(tapemag)
                
            # elseif joint2dprob.ADkindmag=="FWDdiff"
            #     if length(vecmodmag)<20
            #         chuncksize = length(vecmodmag)
            #     else
            #         chuncksize = 20 # 20
            #     end
            #     cfggrdmag = ForwardDiff.GradientConfig(joint2dprob.magmisf,vecmodmag,
            #                                            ForwardDiff.Chunk{chuncksize}())
            #     joint2dprob.autodiffstuffmag[] = cfggrdmag

            # else
            #     joint2dprob.autodiffstuffmag[] = nothing
                
            # end

            ############### Grav ################

            joint2dprob.autodiffstuffgrav[] = MG2D.precalcADstuffgrav(joint2dprob.gravmisf,joint2dprob.ADkindgrav,
                                                                            vecmodgrav)

            # if joint2dprob.ADkindgrav=="REVdiffTAPE"
            #     joint2dprob.autodiffstuffgrav[] = ReverseDiff.GradientTape(joint2dprob.gravmisf,
            #                                                                vecmodgrav)
                
            # elseif joint2dprob.ADkindgrav=="REVdiffTAPEcomp"
            #     # compiled tape
            #     tapegrav = ReverseDiff.GradientTape(joint2dprob.gravmisf,vecmodgrav)
            #     joint2dprob.autodiffstuffgrav[] = ReverseDiff.compile(tapegrav)
                
            # elseif joint2dprob.ADkindgrav=="FWDdiff"
                
            #     if length(vecmodgrav)<20
            #         chuncksize = length(vecmodgrav)
            #     else
            #         chuncksize = 20 # 20
            #     end
                
            #     cfggrdgrav = ForwardDiff.GradientConfig(joint2dprob.gravmisf,vecmodgrav,
            #                                             ForwardDiff.Chunk{chuncksize}())
            #     joint2dprob.autodiffstuffgrav[] = cfggrdgrav

            # else
            #     joint2dprob.autodiffstuffgrav[] = nothing
                
            # end
            
            ## set first time to false!
            joint2dprob.firsttime_grad[]=false
        end
        
        # compute gradient
        vecgrad = ∇misfjoint(joint2dprob,vecmodmag,vecmodgrav)
        
        # return flattened gradient
        return vecgrad

    else
        error("joint2dprob::Joint2DpolyProb(): Wrong argument 'kind': $kind...")
    end
end

#################################################################

"""
 Perform tests on polygons in case of Joint Potential Field HMC inversion.
"""
function (joint2dprob::Joint2DpolyProb)(mcur::Vector{Float64},mnew::Vector{Float64},pnew::Vector{Float64},
                                        lowcon::Union{Vector{Float64},Nothing},upcon::Union{Vector{Float64},Nothing},
                                        LcholmassM::Union{LowerTriangular{Float64,Matrix{Float64}},Diagonal{Float64}},
                                        dowhat::String)
    if dowhat=="performchecks"
        
        # perform checks
        allgood = hmcpolycheckjoint!(joint2dprob,mcur,mnew,pnew,lowcon,upcon,LcholmassM) 
        
    end
    return allgood
end

##############################################################

function (joint2dprob::Joint2DpolyProb)(dowhat::String ;
                                        newbodyindices::Union{Nothing,Vector{<:Vector{<:Integer}}}=nothing)

    if dowhat=="newinds"
        if newbodyindices==nothing
            joint2dprob.orgbodyindices .= joint2dprob.magmisf.bodyindices
        else
            joint2dprob.orgbodyindices .= newbodyindices
        end
        
    elseif dowhat=="backinds"
        joint2dprob.magmisf.bodyindices  .= joint2dprob.orgbodyindices
        joint2dprob.gravmisf.bodyindices .= joint2dprob.orgbodyindices

    elseif dowhat=="getbodyindices"
        @assert joint2dprob.magmisf.bodyindices==joint2dprob.gravmisf.bodyindices
        # COPY the bodyindices otherwise they'll get modified!
        return copy(joint2dprob.magmisf.bodyindices) 
    end
    
    return 
end

######################################################

function (joint2dprob::Joint2DpolyProb)(fid_h5::HDF5.File,dowhat::String)

    if dowhat=="readh5bodyindices"

        nbo = length(joint2dprob.orgbodyindices)
        for i=1:nbo
            dset = "curbodyindices_$i"
            joint2dprob.orgbodyindices[i] = fid_h5[dset][:,end]
        end

    end
    return 
end

##############################################################

function (joint2dprob::Joint2DpolyProb)(fid_h5::HDF5.File,it::Integer,sit::Integer,dowhat::String)

    if dowhat=="savebodyindices"

        nbo = length(joint2dprob.orgbodyindices)
        if it==1
            for i=1:nbo
                nind = length(joint2dprob.orgbodyindices[i])
                bi_h5 = HDF5.create_dataset(fid_h5,"curbodyindices_$i", Int64, ((nind,1), (nind,-1)),chunk=(nind,1))
                HDF5.h5d_close(bi_h5)
            end
        end
    
        # models
        for i=1:nbo
            
            nind = length(joint2dprob.orgbodyindices[i])
            bi_h5 = HDF5.open_dataset(fid_h5,"curbodyindices_$i")
            HDF5.set_extent_dims(bi_h5,(nind,sit))
            bi_h5[:,sit] = joint2dprob.orgbodyindices[i][:]
            HDF5.h5d_close(bi_h5)

        end

    end
    return 
end

#################################################################

function splitmaggrav(joint2dprob::Joint2DpolyProb,vecmodpar::Vector)

    vecmodmag,vecmodgrav = MG2D.splitmaggrav(joint2dprob.magmisf,
                                             joint2dprob.gravmisf,vecmodpar)


    # nbi = length(joint2dprob.orgbodyindices)

    # if joint2dprob.magmisf.whichpar==:all && joint2dprob.gravmisf.whichpar==:all

    #     # mag
    #     vecmodmag = copy(vecmodpar[1:end-nbi])
    #     # grav
    #     vecmodgrav = vcat(vecmodpar[1:end-(7*nbi)],vecmodpar[end-nbi+1:end])
            
    # elseif joint2dprob.magmisf.whichpar==:all && joint2dprob.gravmisf.whichpar==:vertices

    #     # mag
    #     vecmodmag = copy(vecmodpar[1:end])
    #     # grav
    #     vecmodgrav = copy(vecmodpar[1:end-(6*nbi)])

    # elseif joint2dprob.magmisf.whichpar==:vertices && joint2dprob.gravmisf.whichpar==:all

    #     # mag
    #     vecmodmag = copy(vecmodpar[1:end-nbi])
    #     # grav
    #     vecmodgrav = copy(vecmodpar[1:end])

    # elseif joint2dprob.magmisf.whichpar==:magnetization && joint2dprob.gravmisf.whichpar==:density
        
    #     # mag
    #     vecmodmag = copy(vecmodpar[1:end-nbi])
    #     # grav
    #     vecmodgrav = copy(vecmodpar[end-nbi+1:end])
        
    # elseif joint2dprob.magmisf.whichpar==:vertices && joint2dprob.gravmisf.whichpar==:vertices

    #     # mag
    #     vecmodmag = copy(vecmodpar) 
    #     # grav
    #     vecmodgrav = copy(vecmodpar)
        
    # end
    
    return vecmodmag,vecmodgrav 
end

###########################################################################

function ∇misfjoint(joint2dprob::Joint2DpolyProb,vecmodmag::AbstractArray,vecmodgrav::AbstractArray) #vecmodpar::AbstractArray)
    
    #::Union{ReverseDiff.GradientTape,ReverseDiff.CompiledTape,Nothing}
    # split mag and grav vectors of model parameters
    # vecmodmag,vecmodgrav = splitmaggrav(joint2dprob,vecmodpar)
    
    # compute separate gradients for mag & grav
    gradmag = MG2D.calc∇misfmag(joint2dprob.magmisf,vecmodmag,
                              joint2dprob.ADkindmag,joint2dprob.autodiffstuffmag[])
    
    
    gradgrav = MG2D.calc∇misfgrav(joint2dprob.gravmisf,vecmodgrav,
                                joint2dprob.ADkindgrav,joint2dprob.autodiffstuffgrav[])

    
    # merge mag & grav gradients
    #bodyindices = joint2dprob.mymisfmag.bodyindices
    nbi = length(joint2dprob.orgbodyindices)

    if joint2dprob.magmisf.whichpar==:all && joint2dprob.gravmisf.whichpar==:all
        
        gradver = gradmag[1:end-(6*nbi)]+gradgrav[1:end-nbi]
        gradmgn = gradmag[end-(6*nbi)+1:end]
        graddens= gradgrav[end-nbi+1:end]
        grad = vcat(gradver,gradmgn,graddens)
        
    elseif joint2dprob.magmisf.whichpar==:all && joint2dprob.gravmisf.whichpar==:vertices

        gradver = gradmag[1:end-(6*nbi)]+gradgrav[1:end]
        gradmgn = gradmag[end-(6*nbi)+1:end]
        grad = vcat(gradver,gradmgn)
        
    elseif joint2dprob.magmisf.whichpar==:vertices && joint2dprob.gravmisf.whichpar==:all

        gradver = gradmag[1:end]+gradgrav[1:end-nbi]
        graddens = gradgrav[end-nbi+1:end]
        grad = vcat(gradver,graddens)
        
    elseif joint2dprob.magmisf.whichpar==:magnetization && joint2dprob.gravmisf.whichpar==:density
        
        grad = vcat(gradmag,gradgrav)
        
    elseif joint2dprob.magmisf.whichpar==:vertices && joint2dprob.gravmisf.whichpar==:vertices
        
        grad = gradmag+gradgrav
        
    end      
        
    return grad
end

##############################################################

function jointprob2vec(jointprob::Joint2DpolyProb,jointpbod::JointPolygBodies2D)
    jointstruct2vec(jointprob.magmisf.whichpar,jointprob.gravmisf.whichpar,jointpbod)
end

#########################################################################

"""
$(TYPEDSIGNATURES)

Function carrying out all the checks in the case of model parameterization 
as polygonal bodies in the framework of HMC inversion.     
"""
function hmcpolycheckjoint!(jointprob::Joint2DpolyProb,mcur::Vector{Float64},mnew::Vector{Float64},pnew::Vector{Float64},
                            lowcon::Union{Vector{Float64},Nothing},upcon::Union{Vector{Float64},Nothing},
                            LcholmassM::Union{LowerTriangular{Float64,Matrix{Float64}},Diagonal{Float64}} ;
                            maxclk::Integer=3)

    @assert jointprob.magmisf.bodyindices==jointprob.gravmisf.bodyindices

    ##===========================================================
    ## CHECKS on magnetization, densities, vertices etc. in case of polygons
    ##===========================================================
    # No checks if whichpar = magnetization && density and initial checks already performed
    if !jointprob.firstcheck[] && (jointprob.magmisf.whichpar==:magnetization && jointprob.gravmisf.whichpar==:density) 
        return true
    end

    #-------------------------
    ## create some structures 
    qbody = vecmodpar2jointstruct(jointprob.magmisf,jointprob.gravmisf,jointprob.magmisf.bodyindices,mnew)

    ##-------------------------------------------------
    # Initial check intersections of segments and with topography
    #--> Aborting in ase of intersections
    if jointprob.firstcheck[]
        intersect = GeoPoly.checkall(qbody.geom.bo,jointprob.topography)
        if !all(intersect.==:none)
            error("Possible either interesections between polygons or crossings with topography in your starting model. Re-check your parameterization! Aborting!")
        end
        jointprob.firstcheck[]=false
        return nothing
    end
        
    ## fix problems with topography and intersections
    clk = 0

    ##-------------------------------------------------
    ## check ordering of vertices (must be anticlockwise)
    for k=1:qbody.geom.nbo
        clockw = checkanticlockwiseorder(qbody.geom.bo[k])
        
        if !clockw
            jointprob.magmisf.bodyindices[k] = reverse(jointprob.magmisf.bodyindices[k],dims = 1)
            jointprob.gravmisf.bodyindices[k] = reverse(jointprob.gravmisf.bodyindices[k],dims = 1)
        end
    end

    ## new body
    qbody = JointPolygBodies2D(jointprob.magmisf.bodyindices,qbody.geom.allvert,qbody.Jind,qbody.Jrem,qbody.rho,ylatext=jointprob.magmisf.ylatext)

    #----------------------------------------------------
    # check intersections of segments and with topography
    intersect = GeoPoly.checkall(qbody.geom.bo,jointprob.topography)
    
    if jointprob.trytofixpolygons

        ### attempt to fix polygons
        while !all(intersect.==:none)

            # performing all fixing in a random order
            GeoPoly.fixall!(jointprob.topography,qbody.geom.bo,intersect,jointprob.magmisf.bodyindices,lowcon,upcon)
            
            ##-------------------------------------------------
            ## check ordering of vertices (must be anticlockwise)
            for k=1:qbody.geom.nbo
                clockw = checkanticlockwiseorder(qbody.geom.bo[k])
                
                if !clockw
                    jointprob.magmisf.bodyindices[k] = reverse(jointprob.magmisf.bodyindices[k],dims = 1)
                    jointprob.gravmisf.bodyindices[k] = reverse(jointprob.gravmisf.bodyindices[k],dims = 1)
                end
            end
            
            qbody = JointPolygBodies2D(jointprob.magmisf.bodyindices,qbody.geom.allvert,qbody.Jind,qbody.Jrem,qbody.rho,ylatext=jointprob.magmisf.ylatext)
            
            # do intersection checks
            intersect = GeoPoly.checkall(qbody.geom.bo,jointprob.topography)
            
            clk += 1
            
            if clk > maxclk
                # check succeded? nope...
                jointprob("backinds")
                return false
            end    
        end

    else

        ## do NOT attempt to fix polygons        
        if !all(intersect.==:none)
            # check succeded? nope...
            jointprob("backinds")
            return false
        end

    end

    ##-------------------------------------------------
    ## unroll back to vectors
    mnew .= jointprob2vec(jointprob,qbody)
    
    ## Adjust p to the new altered trajectory
    ## pnew = M * grad(K) [ from gradk(K) = M^-1 * p ]
    ##  where grad(K) = mnew - mcur
    ## M is given by L * L^t, the Cholesky decomposition
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
