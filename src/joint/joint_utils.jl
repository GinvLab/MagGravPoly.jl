


using PrettyTables

##############################################################################

"""
Print the structure of model parameters on the screen.
"""
function printjointmodparinfo(mstart,bodyindices,mag_whichpar,grav_whichpar,mlow,mup)

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
    
    if mag_whichpar==:all && grav_whichpar==:all

        ncoo = Integer((nmodpar-(7*nbo))/2)
        vecxz = vcat(["x$i" for i=1:ncoo],["z$i" for i=1:ncoo])
        vecm = vcat(["JiM$i" for i=1:nbo],["JiI$i" for i=1:nbo],["JiD$i" for i=1:nbo],["JrM$i" for i=1:nbo],["JrI$i" for i=1:nbo],["JrD$i" for i=1:nbo])
        vecg = ["ρ$i" for i=1:nbo]
        vecparl[:,1].=vcat(vecxz,vecm,vecg)
        vecparu[:,1].=vcat(vecxz,vecm,vecg)
        
    elseif mag_whichpar==:all && grav_whichpar==:vertices

        ncoo = Integer((nmodpar-(6*nbo))/2)
        vecxz = vcat(["x$i" for i=1:ncoo],["z$i" for i=1:ncoo])
        vecm = vcat(["JiM$i" for i=1:nbo],["JiI$i" for i=1:nbo],["JiD$i" for i=1:nbo],["JrM$i" for i=1:nbo],["JrI$i" for i=1:nbo],["JrD$i" for i=1:nbo])
        vecparl[:,1].=vcat(vecxz,vecm)
        vecparu[:,1].=vcat(vecxz,vecm)
        
    elseif mag_whichpar==:vertices && grav_whichpar==:all

        ncoo = Integer((nmodpar-nbo)/2)
        vecxz = vcat(["x$i" for i=1:ncoo],["z$i" for i=1:ncoo])
        vecg = ["ρ$i" for i=1:nbo]
        vecparl[:,1].=vcat(vecxz,vecg)
        vecparu[:,1].=vcat(vecxz,vecg)
        
    elseif mag_whichpar==:vertices && grav_whichpar==:vertices

        ncoo = Integer(nmodpar/2)
        vecxz = vcat(["x$i" for i=1:ncoo],["z$i" for i=1:ncoo])
        vecparl[:,1].=vecxz
        vecparu[:,1].=vecxz
        
    elseif mag_whichpar==:magnetization && grav_whichpar==:density

        vecm = vcat(["JiM$i" for i=1:nbo],["JiI$i" for i=1:nbo],["JiD$i" for i=1:nbo],["JrM$i" for i=1:nbo],["JrI$i" for i=1:nbo],["JrD$i" for i=1:nbo])
        vecg = ["ρ$i" for i=1:nbo]
        vecparl[:,1].=vcat(vecm,vecg)
        vecparu[:,1].=vcat(vecm,vecg)
        
    end

    vecparl[:,2].=mlow
    vecparu[:,2].=mup
    
    #println(" ---------------------------------------------------")
    printstyled("\n Lower constrains for model parameters:\n",bold=false,color=:yellow)
    pretty_table(vecparl, noheader = true, crop = :horizontal, formatters = ft_round(3))
    #println("\n ---------------------------------------------------")

    printstyled("\n Upper constrains for model parameters:\n",bold=false,color=:yellow)
    pretty_table(vecparu, noheader = true, crop = :horizontal, formatters = ft_round(3))  
    println()#"\n ---------------------------------------------------")
    
    return nothing
end
#----------------------------------

####################################################

"""
Function to create a customized mass matrix.
"""
function calcMjoint(indices::Vector{<:Vector{<:Integer}},modpar::Vector{Float64},sigma::Vector{Vector{Float64}},
                    whichparmag::Symbol,whichpargrav::Symbol;
                    corrlength::Union{Vector{Vector{Float64}},Nothing}=nothing)  
    
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
    
    if (whichparmag==:all && whichpargrav==:all)||
        (whichparmag==:vertices && whichpargrav==:all)||
        (whichparmag==:all && whichpargrav==:vertices)||
        (whichparmag==:vertices && whichpargrav==:vertices)
        
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
        errstr = "calcMjoint(): Number of model parameters does not match with 'whichparmag' and/or 'whichpargrav'. Aborting"
        if (nx-ncoo)/nbo == 0
            try
                a1 = ( length.(sigma) == Integer.(ones(nbo).*2) )
                a2 = (  whichparmag == whichpargrav == :vertices )
            catch
                error(errstr)
            end
            
        elseif (nx-ncoo)/nbo == 1
             try
                 a1 = ( (length.(sigma[1:nbo]) == Integer.(ones(nbo).*2)) && (length.(sigma[nbo+1:2*nbo]) == Integer.(ones(nbo))) )
                 a2 = (whichparmag == :vertices && whichpargrav == :all)
catch
                error(errstr)
             end
            
        elseif (nx-ncoo)/nbo == 6
            try
                a1 = ( (length.(sigma[1:nbo]) == Integer.(ones(nbo).*2)) && (length.(sigma[nbo+1:2*nbo]) == Integer.(ones(nbo).*6)) )
                a2 = (whichparmag == :all && whichpargrav == :vertices)
            catch
                error(errstr)
            end

        elseif (nx-ncoo)/nbo == 7
            try
                a1 = ( (length.(sigma[1:nbo]) == Integer.(ones(nbo).*2)) && (length.(sigma[nbo+1:2*nbo]) == Integer.(ones(nbo).*7)) )
                a2 = ( whichparmag == whichpargrav == :all )
            catch
                error(errstr)
            end

        else
            a1 = false
            a2 = false 
        end      
        
        if a1==false || a2==false
            error(errstr)
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
            
            # If whichparmag and/or whichpargrav == :all
            if (nx-ncoo)/nbo == 1        
                
                l1,l2 = l12(1,ncoo,nbo)
                M[indphys[collect(l1:l2)[i]]] = sigma[i+nbo][end].^2
                
            elseif (nx-ncoo)/nbo == 6
                
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

            elseif (nx-ncoo)/nbo == 7

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
                l1,l2 = l12(7,ncoo,nbo)
                M[indphys[collect(l1:l2)[i]]] = sigma[i+nbo][7].^2

            end
            
        end    
        
        ############################
        
        ## Covarianze legate ai parametri fisici mag e grav
    elseif whichparmag==:magnetization && whichpargrav==:density

        ### Checks on problem
        @assert length(sigma)==nbo
        @assert length.(sigma) == Integer.(ones(nbo).*7)
        @assert length(corrlength)==nothing

        if nx/nbo == 7

            indphys = collect(diagind(M))
            
            for i=1:nbo
                M[indphys[1:nbo]] = sigma[i][1].^2
                M[indphys[nbo+1:2*nbo]] = sigma[i][2].^2
                M[indphys[(2*nbo)+1:3*nbo]] = sigma[i][3].^2
                M[indphys[(3*nbo)+1:4*nbo]] = sigma[i][4].^2
                M[indphys[(4*nbo)+1:5*nbo]] = sigma[i][5].^2
                M[indphys[(5*nbo)+1:6*nbo]] = sigma[i][6].^2
                M[indphys[(6*nbo)+1:7*nbo]] = sigma[i][7].^2 
            end
            
        else
            error("Number of model parameters does not match with 'whichparmag' and/or 'whichpargrav'. Aborting")
            
        end
    else
        error("calcMjoint(): the following combination for 'whichparmag' and 'whichpargrav' are the only possibile: 
                    - whichparmag==:vertices && whichpargrav==:vertices;
                    - whichparmag==:vertices && whichpargrav==:all;
                    - whichparmag==:all && whichpargrav==:all;
                    - whichparmag==:all && whichpargrav==:vertices;
                    - whichparmag==:magnetization && whichpargrav==:density.
                    Aborting!")
        
    end
    
    return M
end


##########################################################################
