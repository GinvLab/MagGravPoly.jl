###############################################

"""
$(TYPEDSIGNATURES)

Compute the DC-Shift of gravity data. There are two possibilities for `type`:
- :auto
- :ref_obs
See the MagGrav2Dpoly manual for details and explanations.
"""
function grav_dcshift!(tgravobs::Vector{<:Real},type::Symbol;id::Union{Nothing,<:Integer}=nothing)

    if type == :auto
        @assert id == nothing
        tgravobs.-=mean(tgravobs)
    elseif type == :ref_obs
        @assert typeof(id) <: Integer
        @assert id > 0 && id <= length(tgravobs)
        tgravobs.-=tgravobs[id]     
    else
        error("The only possibilities for `type` are :auto and :ref_obs. Aborting")     
    end
    
    return
end

###############################################

#############################################################################

"""
Print the structure of model parameters on the screen.
"""
function printgravmodparinfo(mstart,bodyindices,grav_whichpar,mlow,mup)
    
    nbo=length(bodyindices)
    nmodpar=length(mstart)

    if mlow!=nothing && mup!=nothing 
        @assert length(mlow)==length(mup)
        @assert length(mlow)==nmodpar
    else
        @assert mlow==mup
    end
    
    vecparl=Array{Union{Float64,String,Nothing},2}(undef,nmodpar,2)
    vecparu=Array{Union{Float64,String,Nothing},2}(undef,nmodpar,2)
    
    if grav_whichpar==:all

        ncoo = Integer((nmodpar-nbo)/2)
        vecxz = vcat(["x$i" for i=1:ncoo],["z$i" for i=1:ncoo])
        vecg = ["ρ$i" for i=1:nbo]
        vecparl[:,1].=vcat(vecxz,vecg)
        vecparu[:,1].=vcat(vecxz,vecg)
        
    elseif grav_whichpar==:vertices

        ncoo = Integer(nmodpar/2)
        vecxz = vcat(["x$i" for i=1:ncoo],["z$i" for i=1:ncoo])
        vecparl[:,1].=vecxz
        vecparu[:,1].=vecxz
        
    elseif grav_whichpar==:density

        vecg = ["ρ$i" for i=1:nbo]
        vecparl[:,1].=vecg
        vecparu[:,1].=vecg
        
    end

    vecparl[:,2].=mlow
    vecparu[:,2].=mup
    
    println(" ---------------------------------------------------")
    printstyled(" Lower constraints for model parameters:\n",bold=true,color=:light_cyan)
    pretty_table(vecparl, noheader = true, crop = :horizontal, formatters = ft_round(3))

    printstyled(" Upper constraints for model parameters:\n",bold=true,color=:light_cyan)
    pretty_table(vecparu, noheader = true, crop = :horizontal, formatters = ft_round(3))  
    
    return nothing
end
#-----------------------------------

###############################################################################

"""
Function to create a customized Mass Matrix
"""
function calcMgrav(indices::Vector{<:Vector{<:Integer}},modpar::Vector{<:Float64},sigma::Vector{<:Vector{<:Float64}},whichpargrav::Symbol;
                   corrlength::Union{Vector{<:Vector{<:Float64}},Nothing}=nothing)  

    #------------------------------------------
    function cgaussian(dist,corrlength)
        if maximum(dist)==0.0
            return 1.0
        else
            @assert(corrlength>0.0)
            return exp.(-(dist./corrlength).^2)
        end
    end
    #------------------------------------------
    
    nbo = length(indices)
    nx = length(modpar)
    
    ## Creazione matrice grossa con dimensioni = lunghezza vettore modpar
    M = zeros(nx,nx)
    #ind = Integer.(collect(1:nx))
 
    if whichpargrav==:all || whichpargrav==:vertices

        @assert corrlength!=nothing
        @assert length(corrlength)==nbo
        @assert length.(corrlength) == Integer.(ones(nbo).*2)
        
        indlast = Integer(0)
        
        if nbo > 1
            indbo = let
                num = length(indices[1])
                for f=2:nbo
                    nbo2 = collect(1:nbo)
                    deleteat!(nbo2, findall(x->x>=f,nbo2))
                    for h in nbo2
                        match=Int64[]
                        for k=1:length(indices[f])
                            if (indices[h][indices[h].==indices[f][k]]) != []
                                append!(match,indices[h][indices[h].==indices[f][k]])
                            end
                        end
                        num += length(indices[f])-length(match)  
                    end
                    num = num
                end
                num = num
            end
        else
            indbo = length(indices[1])
        end
        
        ncoo = indbo*2
        
        ### Checks on problem
        if (nx-ncoo)/nbo == 0
            @assert length(sigma)==nbo
            @assert length.(sigma) == Integer.(ones(nbo).*2)
            @assert whichpargrav == :vertices 
        elseif (nx-ncoo)/nbo == 1
            @assert length(sigma)==nbo*2
            @assert (length.(sigma[1:nbo]) == Integer.(ones(nbo).*2)) && (length.(sigma[nbo+1:2*nbo]) == Integer.(ones(nbo)))
            @assert whichpargrav == :all
        else
            error("Number of model parameters does not match with 'whichpargrav'. Aborting")
        end      
        
        ## Inserimento covarianze su posizioni vertici (x,z) su specifiche corrlength
        indlast = 0
        indphys = collect(diagind(M))
        
        for i=1:nbo
            
            ########################
            
            match = 0
            
            if i>1 && i<=nbo          
                match = let
                    num = 0
                    nbo2 = collect(1:nbo)
                    deleteat!(nbo2, findall(x->x>=i,nbo2))
                    for h in nbo2
                        match=Int64[] 
                        for k=1:length(indices[h])
                            if (indices[i][indices[i].==indices[h][k]]) != []
                                append!(match,indices[i][indices[i].==indices[h][k]])
                            end
                        end
                        num += length(match)
                    end
                    num = num
                end
                num = num
            end

            ind1 = collect(indlast+1:indlast+length(indices[i])-match)
            ind2 = collect(indlast+indbo+1:indlast+indbo+length(indices[i])-match)

            for j=1:length(ind1)

                dist = sqrt.((modpar[ind1].-modpar[ind1[j]]).^2 + (modpar[ind2].-modpar[ind2[j]]).^2)
                
                M[ind1[j],ind1] .= sigma[i][1]^2 .* cgaussian(dist,corrlength[i][1])
                
                M[ind2[j],ind2] .= sigma[i][2]^2 .* cgaussian(dist,corrlength[i][end])
                
            end
            indlast += length(indices[i])-match

            # If whichpargrav == :all
            if (nx-ncoo)/nbo == 1

                M[indphys[collect(ncoo+1:nx)[i]]] = sigma[i+nbo][end].^2
                
            end
            
        end

        ############################
    
        ## Covarianze legate ai parametri fisici grav
    elseif whichpargrav==:density

        ### Checks on problem
        @assert length(sigma)==nbo
        @assert length.(sigma) == Integer.(ones(nbo))
        @assert corrlength==nothing
        
        if nx/nbo == 1

            indphys = collect(diagind(M))
            
            for i=1:nbo
                M[indphys[collect(1:nbo)[i]]] = sigma[i][end].^2
            end
                        
        else
            error("Number of model parameters does not match with 'whichpargrav'. Aborting")
            
        end
        
    else
        error("calcMgrav(): whichpargrav must be 'vertices', 'density' or 'all'. Aborting!")
        
    end
    
    return M
end

###############################################################################
