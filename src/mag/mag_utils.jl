
#########################################################################

"""
$(TYPEDSIGNATURES)

Vector addition of magnetic (remnant + induced) components.
"""
function magcomp(modJind::Real,Iind::Real,Dind::Real,modJrem::Real,Irem::Real,
                 Drem::Real,C::Real)
    
    ## Induced magnetization components
    Jix = modJind*cos(Iind)*cos(C-Dind)
    Jiy = modJind*cos(Iind)*sin(C-Dind)
    Jiz = modJind*sin(Iind)

    ## Remnant magnetization components
    Jrx = modJrem*cos(Irem)*cos(C-Drem)
    Jry = modJrem*cos(Irem)*sin(C-Drem)
    Jrz = modJrem*sin(Irem)

    ## Vector addition    
    Jtotx = Jix+Jrx
    Jtoty = Jiy+Jry
    Jtotz = Jiz+Jrz
   
    return Jtotx,Jtoty,Jtotz
end

##############################################

"""
$(TYPEDSIGNATURES)

Convert from the field H (A/m) to B (nT).
"""
@inline function convert_H_to_B_nT( H_Am::Real ) 
    ## permeabilita' del vuoto 
    ## muzero = 4.0 * pi * 10.0^-7
    ## B nanoTesla
    ## B_nT = ( muzero * H_Am ) * 10.0^9
    B_nT =  pi * 400.0 * H_Am 
    return B_nT
end

###############################################

"""
$(TYPEDSIGNATURES)

Convert from the field B (nT) to H (A/m).
"""
@inline function convert_B_nT_to_H( B_nT::Real ) 
    H_Am = B_nT / (pi * 400.0)
    return H_Am
end

###############################################

"""
$(TYPEDSIGNATURES)

Create a MagnetizVector of Induced Magnetization through the Earth's Magnetic Field parameters (i.e., module, inclination, declination) and the Magnetic Susceptibility. Units are in SI (Earth's Magnetic Field module is expressed in nT).
"""
function create_magind(mod::Real,incl::Real,decl::Real,susc::AbstractVector{<:Real}) 

    nbody = length(susc)
    modmag = (mod/100.0).*susc
    inclmag = (0.0.*susc).+incl
    declmag = (0.0.*susc).+decl

    magind = MagnetizVector(mod=modmag,Ideg=inclmag,Ddeg=declmag)
    
    return magind
end

###############################################

"""
$(TYPEDSIGNATURES)

Compute the DC-Shift of Total-field Magnetic intensity Anomaly (TMA) data. There are two possibilities for `type`:
- :auto
- :ref_obs
See the Mag2Dpoly manual for details and explanations.
"""
function mag_dcshift!(tmagobs::Vector{<:Real},type::Symbol;id::Union{Nothing,<:Integer}=nothing)

    if type == :auto
        @assert id == nothing
        tmagobs.-=mean(tmagobs)
    elseif type == :ref_obs
        @assert typeof(id) <: Integer
        @assert id > 0 && id <= length(tmagobs)
        tmagobs.-=tmagobs[id]     
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
function printmagmodparinfo(mstart,bodyindices,mag_whichpar,mlow,mup)
    
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
    
    if mag_whichpar==:all

        ncoo = Integer((nmodpar-(6*nbo))/2)
        vecxz = vcat(["x$i" for i=1:ncoo],["z$i" for i=1:ncoo])
        vecm = vcat(["JiM$i" for i=1:nbo],["JiI$i" for i=1:nbo],["JiD$i" for i=1:nbo],["JrM$i" for i=1:nbo],["JrI$i" for i=1:nbo],["JrD$i" for i=1:nbo])
        vecparl[:,1].=vcat(vecxz,vecm)
        vecparu[:,1].=vcat(vecxz,vecm)
        
    elseif mag_whichpar==:vertices

        ncoo = Integer(nmodpar/2)
        vecxz = vcat(["x$i" for i=1:ncoo],["z$i" for i=1:ncoo])
        vecparl[:,1].=vecxz
        vecparu[:,1].=vecxz
        
    elseif mag_whichpar==:magnetization

        vecm = vcat(["JiM$i" for i=1:nbo],["JiI$i" for i=1:nbo],["JiD$i" for i=1:nbo],["JrM$i" for i=1:nbo],["JrI$i" for i=1:nbo],["JrD$i" for i=1:nbo])
        vecparl[:,1].=vecm
        vecparu[:,1].=vecm
        
    end

    vecparl[:,2].=mlow
    vecparu[:,2].=mup
    
    println(" ---------------------------------------------------")
    printstyled("\n Lower constrains for model parameters:\n",bold=true,color=:light_cyan)
    pretty_table(vecparl, noheader = true, crop = :horizontal, formatters = ft_round(3))
    println("\n ---------------------------------------------------")

    printstyled("\n Upper constrains for model parameters:\n",bold=true,color=:light_cyan)
    pretty_table(vecparu, noheader = true, crop = :horizontal, formatters = ft_round(3))  
    println("\n ---------------------------------------------------")
    
    return nothing
end
#-----------------------------------

#########################################################################3

"""
Function to create a customized Mass Matrix
"""
function calcMmag(indices::Vector{<:Vector{<:Integer}},modpar::Vector{<:Float64},sigma::Vector{<:Vector{<:Float64}},whichparmag::Symbol;
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
    
    if whichparmag==:all || whichparmag==:vertices

        ## Controlli sugli input
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
            @assert whichparmag == :vertices
        elseif (nx-ncoo)/nbo == 6
            @assert length(sigma)==nbo*2
            @assert (length.(sigma[1:nbo]) == Integer.(ones(nbo).*2)) && (length.(sigma[nbo+1:2*nbo]) == Integer.(ones(nbo).*6))
            @assert whichparmag == :all 
        else
            error("Number of model parameters does not match with 'whichparmag'. Aborting")
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

            #----------------------------------------------
            l12(n,ncoo,nbo) = (ncoo+1+(n-1)*nbo,ncoo+n*nbo)
            #----------------------------------------------
            
            # If whichparmag == :all
            if (nx-ncoo)/nbo == 6

                l1,l2 = l12(1,ncoo,nbo)
                M[indphys[collect(l1:l2)[i]]] = sigma[i+nbo][1].^2
                l1,l2 = l12(2,ncoo,nbo)
                M[indphys[collect(l1:l2)[i]]] = sigma[i+nbo][2].^2
                l1,l2 = l12(3,ncoo,nbo)
                M[indphys[collect(l1:l2)[i]]] = sigma[i+nbo][3].^2
                l1,l2 = l12(4,ncoo,nbo)
                M[indphys[collect(l1:l2)[i]]] = sigma[i+nbo][4].^2
                l1,l2 = l12(5,ncoo,nbo)
                M[indphys[collect(l1:l2)[i]]] = sigma[i+nbo][5].^2
                l1,l2 = l12(6,ncoo,nbo)
                M[indphys[collect(l1:l2)[i]]] = sigma[i+nbo][6].^2
                
            end

        end
        
        ############################
        
        ## Covarianze legate ai parametri fisici mag
    elseif whichparmag==:magnetization

        ### Checks on problem
        @assert length(sigma)==nbo
        @assert length.(sigma) == Integer.(ones(nbo).*6)
        @assert corrlength==nothing
        
        if nx/nbo == 6
            
            indphys = collect(diagind(M))

            for i=1:nbo
                M[indphys[collect(1:nbo)[i]]] = sigma[i][1].^2
                M[indphys[collect(nbo+1:2*nbo)[i]]] = sigma[i][2].^2
                M[indphys[collect((2*nbo)+1:3*nbo)[i]]] = sigma[i][3].^2
                M[indphys[collect((3*nbo)+1:4*nbo)[i]]] = sigma[i][4].^2
                M[indphys[collect((4*nbo)+1:5*nbo)[i]]] = sigma[i][5].^2
                M[indphys[collect((5*nbo)+1:6*nbo)[i]]] = sigma[i][6].^2
            end
            
        else
            error("Number of model parameters does not match with 'whichparmag'. Aborting")
            
        end
        
    else
        error("calcMmag(): whichparmag must be 'vertices', 'magnetization' or 'all'. Aborting!")
        
    end
    
    return M
end

#################################################################################
