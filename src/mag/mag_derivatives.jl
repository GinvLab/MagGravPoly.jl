
########################################
# structure for misfit calculations
######################################################
# struct Mag2DPolyMisf
#     bodyindices::Vector{Vector{<:Integer}}
#     northxax::Real
#     xzobs::Array{<:Real,2}
#     tmagobs::Vector{<:Real}
#     invcovmat::AbstractMatrix{<:Real}
#     whichpar::Symbol
#     allvert::Union{Nothing,Array{Float64,2}}
#     Jind::Union{Nothing,MagnetizVector}
#     Jrem::Union{Nothing,MagnetizVector}
#     ylatext::Union{Nothing,Vector{<:Real}}

    """
$(TYPEDSIGNATURES)

Julia structure to contain all data required for both misfit and misfit gradient calculation.
The `whichpar` Symbol indicates which paramters the user would like to invert. It should be `:all`, `:vertices` or `:magnetization`.  
"""
struct Mag2DPolyMisf{I<:Integer,F<:Real}
    bodyindices::Vector{Vector{<:I}}
    northxax::F
    xzobs::Array{<:F,2}
    tmagobs::Vector{<:F}
    invcovmat::AbstractMatrix{<:F}
    whichpar::Symbol
    allvert::Union{Nothing,Array{Float64,2}}
    Jind::Union{Nothing,MagnetizVector}
    Jrem::Union{Nothing,MagnetizVector}
    ylatext::Union{Nothing,Vector{<:F}}

    function Mag2DPolyMisf(bodyindices::Vector{<:Vector{<:Integer}},northxax::Real,xzobs::Array{<:Real,2},
                           tmagobs::Vector{<:Real},invcovmat::AbstractMatrix{<:Real},whichpar::Symbol ;
                           allvert=nothing,Jind=nothing,Jrem=nothing,
                           ylatext::Union{Nothing,Vector{<:Real}}=nothing)
        
        # magpbod_copy is a struct to be used as a temporary/mutable thing just to
        #   speed up calculations
        # magpbod_copy = deepcopy(magpbod)
        bodyindices_copy = deepcopy(bodyindices)
        nbo=length(bodyindices)
        
        if whichpar==:all
            allvert=nothing
            Jind=nothing
            Jrem=nothing

        elseif whichpar==:vertices
            allvert=nothing
            @assert Jind!=nothing
            @assert Jrem!=nothing
            @assert length(Jind.mod)==length(Jrem.mod)==nbo
            
        elseif whichpar==:magnetization
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
            @assert length(unidxs)==length(allvert[:,1])

            Jind=nothing
            Jrem=nothing
        else
            error("Mag2DPolyMisf(): 'whichpar' must be ':all', ':vertices' or ':magnetization'. Aborting!")
        end

        return new{typeof(bodyindices_copy[1][1]),typeof(xzobs[1,1])}(bodyindices_copy,northxax,
                                                                      xzobs,tmagobs,invcovmat,
                                                                      whichpar,allvert,Jind,Jrem,ylatext)
        # return new(bodyindices_copy,northxax,xzobs,tmagobs,invcovmat,
        #            whichpar,allvert,Jind,Jrem,ylatext)
    end
end

################################################################

"""
$(TYPEDSIGNATURES)

Function to compute the value of the misfit functional with respect to the model parameters.
"""
function calcmisfmag(modpar::AbstractArray,magmisf::Mag2DPolyMisf)
    misf = magmisf(modpar)
    return misf
end

################################################################

function (magmisf::Mag2DPolyMisf)(modpar::AbstractArray)

    nbo = length(magmisf.bodyindices)
    
    if magmisf.whichpar==:all
        
        ncoo = length(modpar)-(6*nbo)
        nvert = div(ncoo,2)
       
        # set vertices
        allvert = reshape(modpar[1:ncoo],nvert,2)
        
        l12(n,ncoo,nbo) = (ncoo+1+(n-1)*nbo,ncoo+n*nbo)
        # set induced magnetization 
        l1m,l2m = l12(1,ncoo,nbo)
        l1I,l2I = l12(2,ncoo,nbo)
        l1D,l2D = l12(3,ncoo,nbo)
        Jind_mod = modpar[l1m:l2m]
        Jind_Ideg = modpar[l1I:l2I]
        Jind_Ddeg = modpar[l1D:l2D]

        # set remnant magnetization
        l1m,l2m = l12(4,ncoo,nbo)
        l1I,l2I = l12(5,ncoo,nbo)
        l1D,l2D = l12(6,ncoo,nbo)
        Jrem_mod = modpar[l1m:l2m]
        Jrem_Ideg = modpar[l1I:l2I]
        Jrem_Ddeg = modpar[l1D:l2D]
        
    elseif magmisf.whichpar==:vertices

        ncoo = length(modpar) 
        nvert = div(ncoo,2)

        # set vertices
        allvert = reshape(modpar[1:ncoo],nvert,2)

        # set induced magnetization
        Jind_mod = magmisf.Jind.mod
        Jind_Ideg = magmisf.Jind.Ideg
        Jind_Ddeg = magmisf.Jind.Ddeg

        # set remnant magnetization
        Jrem_mod = magmisf.Jrem.mod
        Jrem_Ideg = magmisf.Jrem.Ideg
        Jrem_Ddeg = magmisf.Jrem.Ddeg

    elseif magmisf.whichpar==:magnetization

        # set vertices
        allvert = magmisf.allvert

        # set induced magnetization
        Jind_mod  = modpar[1:nbo]
        Jind_Ideg = modpar[nbo+1:2*nbo]
        Jind_Ddeg = modpar[(2*nbo)+1:3*nbo]

        # set remnant magnetization
        Jrem_mod = modpar[(3*nbo)+1:4*nbo]
        Jrem_Ideg = modpar[(4*nbo)+1:5*nbo]
        Jrem_Ddeg = modpar[(5*nbo)+1:6*nbo]

    end

    #tmag = zeros(eltype(Jinds[1].mod),size(xzobs,1))
    tmagAD = zeros(eltype(modpar),size(magmisf.xzobs,1))

    if magmisf.ylatext!=nothing   
        for i=1:nbo
            curbo = BodySegments2D(magmisf.bodyindices[i],allvert)
            tmagAD .+= tmagpoly2_75D(magmisf.xzobs,Jind_mod[i],Jind_Ideg[i],Jind_Ddeg[i],
                                    Jrem_mod[i],Jrem_Ideg[i],Jrem_Ddeg[i],
                                    magmisf.northxax,curbo,magmisf.ylatext[1],magmisf.ylatext[end]) 
        end
    else    
        for i=1:nbo
            curbo = BodySegments2D(magmisf.bodyindices[i],allvert)
            tmagAD .+= tmagpoly2D(magmisf.xzobs,Jind_mod[i],Jind_Ideg[i],Jind_Ddeg[i],
                                  Jrem_mod[i],Jrem_Ideg[i],Jrem_Ddeg[i],
                                  magmisf.northxax,curbo) 
        end
    end
    
    ##----------------
    # println("misfit")
    # @time begin
    #dcshift = mean(magmisf.tmagobs)-mean(tmagAD)
    #tmagAD.+=dcshift
    #tmagobs=magmisf.tmagobs.-dcshift
    #dif=tmagAD.-tmagobs
    dif = tmagAD.-magmisf.tmagobs
    tmp = magmisf.invcovmat * dif 
    misf = 0.5 .* dot(dif,tmp)
    # end
    
    return misf
end

#######################################################
"""
$(TYPEDSIGNATURES)

Pre-calculate some parameters before computing the gradient of the misfit using automatic differentiation.

"""
function precalcADstuffmag(magmisf::Mag2DPolyMisf,ADkind::String,vecmodpar::AbstractArray)

    if ADkind=="REVdiffTAPE"
        autodiffstuff = ReverseDiff.GradientTape(magmisf,vecmodpar)

    elseif ADkind=="REVdiffTAPEcomp"
        # compiled tape
        tape = ReverseDiff.GradientTape(magmisf,vecmodpar)
        autodiffstuff = ReverseDiff.compile(tape)

    elseif ADkind=="FWDdiff"
        if length(vecmodpar)<20
            chuncksize = length(vecmodpar)
        else
            chuncksize = 20 # 20
        end
        cfggrd = ForwardDiff.GradientConfig(magmisf,vecmodpar,ForwardDiff.Chunk{chuncksize}())
        autodiffstuff = cfggrd

    else
        autodiffstuff = nothing
        
    end

    return autodiffstuff
end


########################################
"""
$(TYPEDSIGNATURES)

Function to compute the gradient of misfit with respect to the model parameters required for HMC inversions. 
The gradient is computed by means of automatic differentiation using one of three different methods. The user must indicate the method by the `ADkind` variable, choosing among `FWDdiff`, `REVdiffTAPE` or `REVdiffTAPEcomp`. 
For an explanation about the automatic differentiation method, the reader is invited to look at the documentation 
relative to the Julia package `ForwardDiff` and `ReverseDiff`.
"""  
function calc∇misfmag(magmisf::Mag2DPolyMisf,modpar::AbstractArray, # ,whichpar::String,
                      ADkind::String,autodiffstuff )    

    if ADkind=="REVdiffTAPE"
        #println("\nRev diff tape")
        grad = similar(modpar)
        ReverseDiff.gradient!(grad,autodiffstuff,modpar)

     elseif ADkind=="REVdiffTAPEcomp"
        #println("\nRev diff tape")
        grad = similar(modpar)
        ReverseDiff.gradient!(grad,autodiffstuff,modpar)

     elseif ADkind=="FWDdiff"
        #println("\nFWD diff")
        # chuncksize = 20 # 20
        # cfggrd = ForwardDiff.GradientConfig(magmisf,modpar,ForwardDiff.Chunk{chuncksize}())
        # allocate output
        grad = similar(modpar)
        ForwardDiff.gradient!(grad,magmisf,modpar,autodiffstuff)

    end

    nangrad = isnan.(grad)
    if any(nangrad)
        error("∇misf(): magmisf error, isnan.(grad) = $nangrad")
    end
    return grad
end

#######################################################################

function jacfindiff(xzobs,northxax,pbody; δx=0.1)    
    # Jacobian by finite difference

    nobs = size(xzobs,1)
    nver = length(pbody.geom.allvert)
    nmpar = nver+6*pbody.geom.nbo
    Jfd = zeros(nobs,nmpar)

    tm0 = tmagpolybodies2D(xzobs,northxax,pbody)

    for o=1:nobs
        for v=1:nver
            pbody1 = deepcopy(pbody)
            pbody1.geom.allvert[v] += δx
            tm1 = tmagpolybodies2D(xzobs[o:o,:],northxax,pbody1)
            Jfd[o,v] = (tm1[1].-tm0[o])./δx
        end

        l=nver+1

        for m=1:6

            for j=1:pbody.geom.nbo

                pbody1 = deepcopy(pbody)

                if m==1
                    pbody1.Jind.mod[j] += δx
                elseif m==2
                    pbody1.Jind.Ideg[j] += δx
                elseif m==3
                    pbody1.Jind.Ddeg[j] += δx
                elseif m==4
                    pbody1.Jrem.mod[j] += δx
                elseif m==5
                    pbody1.Jrem.Ideg[j] += δx
                elseif m==6
                    pbody1.Jrem.Ddeg[j] += δx
                end

                tm1 = tmagpolybodies2D(xzobs[o:o,:],northxax,pbody1)
                Jfd[o,l] = (tm1[1].-tm0[o])./δx
                l+=1
            end

        end

    end
    return Jfd
end

###########################################################################

