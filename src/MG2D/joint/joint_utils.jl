
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
    
    println(" ---------------------------------------------------")
    printstyled(" Lower constraints for model parameters:\n",bold=true,color=:light_cyan)
    pretty_table(vecparl, formatters = [fmt__round(3)])
    #pretty_table(vecparl, noheader = true, crop = :horizontal, formatters = ft_round(3))

    printstyled(" Upper constraints for model parameters:\n",bold=true,color=:light_cyan)
    pretty_table(vecparu, formatters = [fmt__round(3)])
    #pretty_table(vecparu, noheader = true, crop = :horizontal, formatters = ft_round(3))  
    
    return nothing
end
#----------------------------------

####################################################
#=
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

        @show (nx-ncoo)/nbo,nx,ncoo,nbo,indbo

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

=#
##########################################################################

"""
$(TYPEDSIGNATURES)

Function to create a customized mass matrix for joint magnetic–gravity problems.

Constructs a full covariance/mass matrix `M` containing:
- Vertex coordinate covariances (x, z);
- Optional physical parameter variances (magnetization and/or density).

# Arguments
- `indices`        :: Vector of vertex index lists per body
- `modpar`         :: Vector of model parameters [x₁..xₙ, z₁..zₙ, (phys params...)]
- `sigma`          :: Vector of per-body standard deviations
- `whichparmag`    :: Symbol — :vertices, :all, :magnetization
- `whichpargrav`   :: Symbol — :vertices, :all, :density
- `corrlength`     :: Vector of per-body correlation lengths (or `nothing`)

# Returns
A symmetric matrix `M` of size `(nx, nx)` representing the joint covariance structure and its inverse.
"""
function calcMjoint(indices::Vector{<:Vector{<:Integer}},
                    modpar::Vector{Float64},
                    sigma::Vector{Vector{Float64}},
                    whichparmag::Symbol,
                    whichpargrav::Symbol;
                    corrlength::Union{Vector{Vector{Float64}},Nothing}=nothing)

    # -------------------------------
    # Internal helper: Gaussian correlation
    # -------------------------------
    cgaussian(dist, corrlength) =  exp.(-(dist ./ corrlength) .^ 2)

    # -------------------------------
    # Internal helper: Safe way to calc. inv
    # -------------------------------
    function invposdefmat(A)
        if !isposdef(A)
            A .+=  I(size(A,1)).*10^-10
        end
        L = cholesky(A).L  ## C = L L'
        L_inv = inv(L)
        A_inv = L_inv'*L_inv ##  (AB)^-1 = B^-1 A^-1
        return A_inv
    end

    # -------------------------------
    # Step 1: Initialize quantities
    # -------------------------------
    nbo = length(indices)
    nx  = length(modpar)
    M   = zeros(nx, nx)
    Minv = zeros(nx, nx)

    # Build compact vertex mapping (unique vertex indexing)
    all_idx = collect(Iterators.flatten(indices))
    unique_list = unique(all_idx)
    indbo = length(unique_list)
    ncoo  = indbo * 2

    if nx < ncoo
        error("calcMjoint(): modpar too short to contain vertex coordinates.")
    end

    xcoords = modpar[1:indbo]
    zcoords = modpar[indbo+1:2*indbo]
    pos_of_old = Dict(old => i for (i, old) in enumerate(unique_list))

    # -------------------------------
    # Step 2: Case A — Vertex covariances (vertices or all)
    # -------------------------------
    if (whichparmag in (:vertices, :all) || whichparmag == :all) &&
       (whichpargrav in (:vertices, :all) || whichpargrav == :all)

        @assert corrlength !== nothing "calcMjoint(): corrlength required for :vertices or :all"
        @assert length(corrlength) == nbo "calcMjoint(): corrlength must have one entry per body"

        # Check shapes
        for i in 1:nbo
            @assert length(corrlength[i]) >= 2 "corrlength[$i] must have [x,z]"
            @assert length(sigma[i]) >= 2 "sigma[$i] must have [σx, σz]"
        end

        # Compact vertex indices per body
        body_pos = [ [pos_of_old[v] for v in inds] for inds in indices ]

        # -------------------------------
        # Step 3: Build vertex covariance submatrices
        # -------------------------------
        for i in 1:nbo
            poslist = body_pos[i]
            ni = length(poslist)
            if ni == 0; continue; end

            xi, zi = xcoords[poslist], zcoords[poslist]
            Xi, Zi = reshape(xi, (ni,1)), reshape(zi, (ni,1))
            distmat = sqrt.((Xi .- Xi').^2 .+ (Zi .- Zi').^2)

            ind_x = poslist
            ind_z = poslist .+ indbo
            
            σx, σz = sigma[i][1], sigma[i][2]
            clx, clz = corrlength[i][1], corrlength[i][end]

            M[ind_x, ind_x] .= σx^2 .* cgaussian(distmat, clx) #clx
            M[ind_z, ind_z] .= σz^2 .* cgaussian(distmat, clz) #clz

            # Calculating the inverse...
            if whichparmag == :vertices && whichpargrav == :vertices
                Minv[ind_x, ind_x] .= invposdefmat(M[ind_x, ind_x])
                Minv[ind_z, ind_z] .= invposdefmat(M[ind_z, ind_z])
            end       
        end
        
        # -------------------------------
        # Step 4: Add physical parameter variances (magnetic + gravity)
        # -------------------------------
        nbphys_total = nx - ncoo
        nbphys_per_body = nbphys_total ÷ max(nbo, 1)

        # Cases:
        #   (nx - ncoo)/nbo == 0  → only vertex coords
        #   == 1                 → gravity density
        #   == 6                 → magnetization
        #   == 7                 → both

        if nbphys_per_body == 0
            # Only vertex coordinates
            nothing

        elseif nbphys_per_body == 1
            # Magnetic vertices + gravity densities
            @assert whichparmag == :vertices && whichpargrav == :all ||
                    whichparmag == :all && whichpargrav == :vertices ||
                    whichparmag == :all && whichpargrav == :all "calcMjoint(): inconsistent param combination for nbphys=1"
            @assert length(sigma) >= 2*nbo "calcMjoint(): sigma must contain entries for physical params"
            for i in 1:nbo
                pos = 2*indbo + i
                M[pos, pos] = sigma[nbo + i][end]^2
            end

        elseif nbphys_per_body == 6
            # Gravity vertices + magnetization
            @assert whichparmag == :all && whichpargrav == :vertices ||
                    whichparmag == :all && whichpargrav == :all "calcMjoint(): inconsistent param combination for nbphys=6"
            @assert length(sigma) >= 2*nbo "calcMjoint(): sigma must include magnetization entries"
            for p in 1:6
                for i in 1:nbo
                    pos = 2*indbo + (p - 1)*nbo + i
                    M[pos, pos] = sigma[nbo + i][p]^2
                end
            end

        elseif nbphys_per_body == 7
            # Full joint model: magnetization + density
            @assert whichparmag == :all && whichpargrav == :all "calcMjoint(): inconsistent param combination for nbphys=7"
            @assert length(sigma) >= 2*nbo "calcMjoint(): sigma must include full magnetization+density entries"
            for p in 1:7
                for i in 1:nbo
                    pos = 2*indbo + (p - 1)*nbo + i
                    M[pos, pos] = sigma[nbo + i][p]^2
                end
            end

        else
            error("calcMjoint(): invalid model layout — (nx - ncoo)/nbo = $nbphys_per_body not supported.")
        end

        # Calculating the inverse...
        if whichparmag != :vertices && whichpargrav != :vertices
            Minv = invposdefmat(M)
        end
        
    # -------------------------------
    # Step 5: Case B — Magnetization + Density only
    # -------------------------------
    elseif whichparmag == :magnetization && whichpargrav == :density
        @assert corrlength === nothing "calcMjoint(): corrlength must be nothing for pure physical param case"
        @assert length(sigma) == nbo "calcMjoint(): sigma must have nbo entries"
        for i in 1:nbo
            @assert length(sigma[i]) == 7 "calcMjoint(): sigma[$i] must have 7 entries (magnetization + density)"
        end

        if nx % nbo != 0
            error("calcMjoint(): nx not divisible by nbo in magnetization/density case.")
        end
        nbphys_per_body = nx ÷ nbo
        @assert nbphys_per_body == 7 "calcMjoint(): expected 7 parameters per body."

        for p in 1:7
            for i in 1:nbo
                pos = (p - 1)*nbo + i
                M[pos, pos] = sigma[i][p]^2
            end
        end
        #Calculating the inverse...
        Minv = invposdefmat(M)

    else
        error("""calcMjoint(): invalid combination of whichparmag/whichpargrav.
Allowed combinations:
  - (:vertices, :vertices)
  - (:vertices, :all)
  - (:all, :vertices)
  - (:all, :all)
  - (:magnetization, :density)
Aborting!""")
    end

    # Last check M positive definite
    if !isposdef(M)
        M .+=  I(size(M,1)).*10^-10
    end

    #Last check Minv positive definite
    if !isposdef(Minv)
        Minv .+=  I(size(Minv,1)).*10^-10
    end
        
    return M,Minv
end


############################################
             
       
