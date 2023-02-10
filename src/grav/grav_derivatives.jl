
########################################
##   structure for misfit calculations
###################################################################

# struct Grav2DPolyMisf{I<:Integer,F<:Real}
#     bodyindices::Vector{Vector{I}}
#     xzobs::Array{F,2}
#     tgravobs::Vector{F}
#     invcovmat::AbstractMatrix{F}
#     whichpar::Symbol
#     allvert::Union{Nothing,Array{F,2}}
#     rho::Union{Nothing,Vector{F}}
#     ylatext::Union{Nothing,Vector{F}}

# struct Grav2DPolyMisf 
#     bodyindices::Vector{Vector{<:Integer}}
#     xzobs::Array{<:Real,2}
#     tgravobs::Vector{<:Real}
#     invcovmat::AbstractMatrix{<:Real}
#     whichpar::Symbol
#     allvert::Union{Nothing,Array{Float64,2}}
#     rho::Union{Nothing,Vector{Float64}}
#     ylatext::Union{Nothing,Vector{Float64}}

"""
$(TYPEDSIGNATURES)

$(TYPEDFIELDS)

Structure containing all data required for both misfit and gradient calculations.
The `whichpar` symbol indicates which paramters the user would like to invert. It should be `:all`, `:vertices` or `:density`.  
"""
struct Grav2DPolyMisf{I<:Integer,F<:Real}
    bodyindices::Vector{Vector{I}}
    xzobs::Array{F,2}
    tgravobs::Vector{F}
    invcovmat::AbstractMatrix{F}
    whichpar::Symbol
    allvert::Union{Nothing,Array{F,2}}
    rho::Union{Nothing,Vector{F}}
    ylatext::Union{Nothing,Vector{F}}

    function Grav2DPolyMisf(bodyindices::Vector{<:Vector{<:Integer}},xzobs::Array{<:Real,2},
                            tgravobs::Vector{<:Real},invcovmat::AbstractMatrix{<:Real},whichpar::Symbol ;
                            allvert=nothing,rho=nothing,ylatext::Union{Nothing,Vector{<:Real}}=nothing)
        
        # magpbod_copy is a struct to be used as a temporary/mutable thing just to
        #   speed up calculations
        # magpbod_copy = deepcopy(magpbod)
        bodyindices_copy = deepcopy(bodyindices)
        nbo=length(bodyindices)
        
        if whichpar==:all
            allvert=nothing
            rho=nothing

        elseif whichpar==:vertices
            allvert=nothing
            @assert rho!=nothing
            @assert length(rho)==nbo
            
        elseif whichpar==:density
            @assert allvert!=nothing
            # if nbo > 1
            #     indbo = let
            #         num = length(bodyindices[1])
            #         for i=2:nbo
            #             nbo2 = collect(1:nbo)
            #             deleteat!(nbo2, findall(x->x>=i,nbo2))
            #             for j in nbo2
            #                 match=Int64[]
            #                 for k=1:length(bodyindices[i])
            #                     if (bodyindices[j][bodyindices[j].==bodyindices[i][k]]) != []
            #                         append!(match,bodyindices[j][bodyindices[j].==bodyindices[i][k]])
            #                     end
            #                 end
            #                 num += length(bodyindices[i])-length(match)  
            #             end
            #             num = num
            #         end
            #         num = num
            #     end
            # else
            #     indbo = length(bodyindices[1])
            # end
            # @assert indbo==length(allvert[:,1])
            
            # Get all the unique body indices
            allind = collect(Iterators.flatten(bodyindices))
            unidxs = unique(allind)

            rho=nothing
        else
            error("Grav2DPolyMisf(): 'whichpar' must be ':all', ':vertices' or ':density'. Aborting!")
        end

        # return new(bodyindices_copy,xzobs,tgravobs,invcovmat,
        #            whichpar,allvert,rho,ylatext)
        return new{typeof(bodyindices_copy[1][1]),typeof(xzobs[1,1])}(bodyindices_copy,
                                                                      xzobs,tgravobs,invcovmat,
                                                                      whichpar,allvert,rho,ylatext)
    end
end

#################################################

"""
$(TYPEDSIGNATURES)

Function to compute the value of the misfit functional with respect to the model parameters.
"""
function calcmisfgrav(modpar::AbstractArray,gravmisf::Grav2DPolyMisf)
    misf = gravmisf(modpar)
    return misf
end

###############################################################
 
function (gravmisf::Grav2DPolyMisf)(modpar::AbstractArray)
    
    nbo = length(gravmisf.bodyindices)
    
    if gravmisf.whichpar==:all

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
        rho = gravmisf.rho

    elseif gravmisf.whichpar==:density

        allvert = gravmisf.allvert
        rho  = modpar[1:end]

    end
    
    #tmag = zeros(eltype(Jinds[1].mod),size(xzobs,1))
    tgravAD = zeros(eltype(modpar),size(gravmisf.xzobs,1))

    if gravmisf.ylatext!=nothing
        for i=1:nbo
            curbo = BodySegments2D(gravmisf.bodyindices[i],allvert)
            tgravAD .+= tgravpoly2_75D(gravmisf.xzobs,rho[i],curbo,gravmisf.ylatext[1],gravmisf.ylatext[2]) 
        end
    else
        for i=1:nbo
            curbo = BodySegments2D(gravmisf.bodyindices[i],allvert)
            tgravAD .+= tgravpoly2D(gravmisf.xzobs,rho[i],curbo) 
        end
    end
    
    ##----------------
    # println("misfit")
    # @time begin
    #dcshift = mean(gravmisf.tgravobs)-mean(tgravAD)
    #tgravAD.+=dcshift
    #tgravobs=gravmisf.tgravobs.-dcshift 
    #dif=tgravAD.-tgravobs
    dif = tgravAD.-gravmisf.tgravobs
    tmp = gravmisf.invcovmat * dif
    misf = 0.5 .* dot(dif,tmp)
    # end

    return misf
end

#######################################################
"""
$(TYPEDSIGNATURES)

Pre-calculate some parameters before computing the gradient of the misfit using automatic differentiation.

"""
function precalcADstuffgrav(gravmisf::Grav2DPolyMisf,ADkind::String,vecmodpar::AbstractArray)

    if ADkind=="REVdiffTAPE"
        autodiffstuff = ReverseDiff.GradientTape(gravmisf,vecmodpar)

    elseif ADkind=="REVdiffTAPEcomp"
        # compiled tape
        tape = ReverseDiff.GradientTape(gravmisf,vecmodpar)
        autodiffstuff = ReverseDiff.compile(tape)

    elseif ADkind=="FWDdiff"
        if length(vecmodpar)<20
            chuncksize = length(vecmodpar)
        else
            chuncksize = 20 # 20
        end
        cfggrd = ForwardDiff.GradientConfig(gravmisf,vecmodpar,ForwardDiff.Chunk{chuncksize}())
        autodiffstuff = cfggrd

    else
        autodiffstuff = nothing
        
    end

    return autodiffstuff
end

#######################################################
"""
$(TYPEDSIGNATURES)

Function to compute the gradient of the misfit functional with respect to the model parameters as required for HMC inversions. 
The gradient is computed by means of automatic differentiation using one of three different methods. The user must indicate the method by the `ADkind` variable, choosing among `FWDdiff`, `REVdiffTAPE` or `REVdiffTAPEcomp`. 
For an explanation about the automatic differentiation method, the reader is invited to look at the documentation 
relative to the Julia package `ForwardDiff` and `ReverseDiff`.  
"""  
function calc∇misfgrav(gravmisf::Grav2DPolyMisf,modpar::AbstractArray,#whichpar::String,
                       ADkind::String, autodiffstuff )

    #::Union{ReverseDiff.GradientTape,ReverseDiff.CompiledTape,Nothing}

    if ADkind=="REVdiffTAPE"
        #println("\nRev diff tape")
        grad = similar(modpar)
        ReverseDiff.gradient!(grad,autodiffstuff,modpar)

     elseif ADkind=="REVdiffTAPEcomp"
        #println("\nRev diff tape")
        grad = similar(modpar)
        ReverseDiff.gradient!(grad,autodiffstuff,modpar)

     elseif ADkind=="FWDdiff"
        # allocate output
        grad = similar(modpar)
        ForwardDiff.gradient!(grad,gravmisf,modpar,autodiffstuff)

    end

    nangrad = isnan.(grad)
    if any(nangrad)
        error("∇misf(): gravmisf error, isnan.(grad) = $nangrad")
    end
    return grad
end


########################################
