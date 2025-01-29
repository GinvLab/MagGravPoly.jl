
#########################################

"""
$(TYPEDSIGNATURES)

Function fixing polygons vertices in order to avoid intersection and topography crossing, changing sign to the corresponding elements in the momentum array. This algorithm is in agreement with the mathematical properties of the Hamiltonian Dynamics.  
"""
function fixall!(topo::TopoEdges,qbo::Vector{BodySegments2D},inters::Vector{Symbol},
                 indices::Union{Vector{Vector{<:Integer}},Nothing}=nothing,
                 lowcon::Union{Vector{Float64},Nothing}=nothing,upcon::Union{Vector{Float64},Nothing}=nothing)
    

    @assert typeof(lowcon)==typeof(upcon)
    nchk = findall(x->x!=:none,inters)
    
    #########################################
    
    if lowcon!=nothing && upcon!=nothing
        @assert indices!=nothing
        nbo = length(indices)
        if nbo > 1
            indbo = let
                num = length(indices[1])
                for i=2:nbo
                    nbo2 = collect(1:nbo)
                    deleteat!(nbo2, findall(x->x>=i,nbo2))
                    for j in nbo2
                        match=Int64[]
                        for k=1:length(indices[i])
                            if (indices[j][indices[j].==indices[i][k]]) != []
                                append!(match,indices[j][indices[j].==indices[i][k]])
                            end
                        end
                        num += length(indices[i])-length(match)  
                    end
                    num = num
                end
                num = num
            end
        else
            indbo = length(indices[1])
        end
        
        indver = indbo*2
        
        blow = Vector{BodySegments2D}(undef,nbo)
        bup = Vector{BodySegments2D}(undef,nbo)
        
        for i=1:nbo
            blow[i] = BodySegments2D(indices[i],reshape(lowcon[1:indver],:,2))
            bup[i] = BodySegments2D(indices[i],reshape(upcon[1:indver],:,2))
        end
        
    else
        blow = nothing
        bup = nothing
    end
    
    ######################################
    
    if !isempty(nchk)
        shuffle!(nchk)
        
        for i in nchk
            if i == 1
                verpolyshift!(qbo,blow,bup)
            elseif i ==2
                verpolyallshift!(qbo,blow,bup)
            elseif i == 3
                vertoposhift!(topo,qbo,blow,bup)
            end      
        end
    end
    return nothing
end

####################################################################################################

"""
$(TYPEDSIGNATURES)

Check if there is any self-intersection considering all polygons beyond to topography crossing, returning three string values.
"""
function checkall(bos::Vector{BodySegments2D},topo::TopoEdges)

    ## Check for any self-intersection considering all polygons
    sign1,sign2 = checkpoly(bos)
    
    ## intersection between polygons and topography
    cross = checktopo(topo,bos)
    if cross
        return [sign1,sign2,:crosstopo]
    end
    return [sign1,sign2,:none]
end

#####################################################################################################

"""
$(TYPEDSIGNATURES)

Check if there is any self-intersection considering all polygons, returning two string values.
"""
function checkpoly(bos::Vector{BodySegments2D})

    sign1 = :none
    sign2 = :none

    ##-------------------
    ## self-intersection for each single polygon
    for i=1:length(bos)
        inter = selfintersectpoly(bos[i])
        if inter
            sign1 = :selfinters
        end
    end

    ##-------------------
    ## intersection between pairs of polygons
    inter = intersectpairpoly(bos)
    if inter
        sign2 = :pairinters
    end
        
    return sign1,sign2 
end

##################################################################################################

"""
$(TYPEDSIGNATURES)

Check if there is any self-intersection for a single polygon, returning a Boolean value (i.e. true or false).
"""
function selfintersectpoly(bo::BodySegments2D)

    nedge = length(bo.ver1[:,1])

    # No check if triangle
    if nedge==3
        return false
    end
    
    polyarray = Array{Array{Float64,1},1}(undef,(nedge+1))   
    
    for n=1:nedge
        polyarray[n] = bo.ver1[n,:]
        if n==nedge
            polyarray[end] = bo.ver1[1,:]
        end
    end

    # define polygon for LibGEOS
    poly = LibGEOS.Polygon([polyarray])
    
    if !LibGEOS.isValid(poly)
        return true
    end
    
    return false
end


#####################################################################################################

"""
$(TYPEDSIGNATURES)

Function to check if there is intersection beetween sides of a polygon, trying to shift polygon vertices in order to avoid intersection. This algorithm shifts the polygon vertices in agreement with the mathematical properties of the Hamiltonian Dynamics.
"""
function verpolyshift!(borg::Vector{BodySegments2D},
                       blow::Union{Vector{BodySegments2D},Nothing}=nothing,
                       bup::Union{Vector{BodySegments2D},Nothing}=nothing)  

    @assert typeof(blow)==typeof(bup)
    #idpolyorg = [Int64[] for a in 1:length(borg)]
    fix = true
    countfix = 0
    sign1 = rand(2:5)*length(borg)
    sign2 = rand(2:5)*length(borg)

    while fix

        if countfix == sign1
            break
        end

        countfix += 1
        #idpoly = deepcopy(idpolyorg)
        bo = deepcopy(borg)
        nbo = shuffle(collect(1:length(bo)))
        fix = false
        int = true        
        countint = 0
        
        while int
            
            if countint == sign2
                fix = true
                break
            end
            
            countint += 1
            int = false
            
            for i in nbo
                num_vert = length(bo[i].ver1[:,1])
                
                for ind_v1=1:(num_vert-2)

                    # No check if triangle
                    if num_vert==3
                        break
                    end
                    
                    p1= [bo[i].ver1[ind_v1,1], bo[i].ver1[ind_v1,2]]
                    q1= [bo[i].ver1[ind_v1+1,1], bo[i].ver1[ind_v1+1,2]]
                    
                    for ind_v2=(ind_v1+2):(num_vert) #taking into account last segment vertice -> first vertice

                        p2= [bo[i].ver1[ind_v2,1], bo[i].ver1[ind_v2,2]]
                        
                        if ind_v1>1
                                                        
                            if ind_v2==num_vert
                                q2= [bo[i].ver1[1,1], bo[i].ver1[1,2]]
                                p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])
                                
                                if p0 != nothing
                                    
                                    int = true
                                    dist1 = sqrt(((p2[1]-q2[1])^2)+((p2[end]-q2[end])^2))
                                    dist2 = sqrt(((p0[1]-p2[1])^2)+((p0[end]-p2[end])^2))
                                    ratio = dist2/dist1
                                    Δx = p0[1]-q2[1]
                                    Δz = p0[end]-q2[end]
                                    xtmp = p0[1]-(Δx*ratio)
                                    ztmp = p0[end]-(Δz*ratio)

                                    if blow != nothing && bup != nothing

                                        if (blow[i].ver1[ind_v2,1] <= xtmp <= bup[i].ver1[ind_v2,1])&&
                                            (blow[i].ver1[ind_v2,2] <= ztmp <= bup[i].ver1[ind_v2,2])
                                           
                                            bo[i].ver1[ind_v2,1] = xtmp
                                            bo[i].ver1[ind_v2,2] = ztmp

                                            #--------------
                                            p2= [bo[i].ver1[ind_v2,1], bo[i].ver1[ind_v2,2]]
                                            p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])
                                            
                                            if p0 == nothing
                                                int = false
                                            end
                                            #--------------
                                        end

                                    else
                                        bo[i].ver1[ind_v2,1] = xtmp
                                        bo[i].ver1[ind_v2,2] = ztmp
                                        
                                        #--------------
                                        p2= [bo[i].ver1[ind_v2,1], bo[i].ver1[ind_v2,2]]
                                        p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])
                                        
                                        if p0 == nothing
                                            int = false
                                        end
                                        #--------------
                                    end
                                end
                                
                            else
                                q2= [bo[i].ver1[ind_v2+1,1], bo[i].ver1[ind_v2+1,2]]
                                p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])
                                
                                if p0 != nothing
                                    
                                    int = true
                                    dist1 = sqrt(((p2[1]-q2[1])^2)+((p2[end]-q2[end])^2))
                                    dist2 = sqrt(((p0[1]-p2[1])^2)+((p0[end]-p2[end])^2))
                                    ratio = dist2/dist1
                                    Δx = p0[1]-q2[1]
                                    Δz = p0[end]-q2[end]
                                    xtmp = p0[1] - (Δx*ratio)
                                    ztmp = p0[end] - (Δz*ratio)
                                    
                                    if blow != nothing && bup != nothing
                                        
                                        if (blow[i].ver1[ind_v2,1] <= xtmp <= bup[i].ver1[ind_v2,1])&&
                                            (blow[i].ver1[ind_v2,2] <= ztmp <= bup[i].ver1[ind_v2,2])
                                            
                                            bo[i].ver1[ind_v2,1] = xtmp
                                            bo[i].ver1[ind_v2,2] = ztmp

                                            #--------------
                                            p2= [bo[i].ver1[ind_v2,1], bo[i].ver1[ind_v2,2]]
                                            p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])
                                            
                                            if p0 == nothing
                                                int = false
                                            end
                                            #--------------
                                        end

                                    else
                                        bo[i].ver1[ind_v2,1] = xtmp
                                        bo[i].ver1[ind_v2,2] = ztmp

                                        #--------------
                                        p2= [bo[i].ver1[ind_v2,1], bo[i].ver1[ind_v2,2]]
                                        p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])
                                        
                                        if p0 == nothing
                                            int = false
                                        end
                                        #--------------
                                    end       
                                end
                            end
                        else
                            
                            if ind_v2==num_vert
                                break
                            else
                                q2= [bo[i].ver1[ind_v2+1,1], bo[i].ver1[ind_v2+1,2]]
                                p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])
                                
                                if p0 != nothing
                                    
                                    int = true
                                    dist1 = sqrt(((p2[1]-q2[1])^2)+((p2[end]-q2[end])^2))
                                    dist2 = sqrt(((p0[1]-p2[1])^2)+((p0[end]-p2[end])^2))
                                    ratio = dist2/dist1
                                    Δx = p0[1]-q2[1]
                                    Δz = p0[end]-q2[end]
                                    xtmp = p0[1] - (Δx*ratio)
                                    ztmp = p0[end] - (Δz*ratio)
                                    
                                    if blow != nothing && bup != nothing
                                        
                                        if (blow[i].ver1[ind_v2,1] <= xtmp <= bup[i].ver1[ind_v2,1])&&
                                            (blow[i].ver1[ind_v2,2] <= ztmp <= bup[i].ver1[ind_v2,2])
                                            
                                            bo[i].ver1[ind_v2,1] = xtmp
                                            bo[i].ver1[ind_v2,2] = ztmp

                                            #--------------
                                            p2= [bo[i].ver1[ind_v2,1], bo[i].ver1[ind_v2,2]]
                                            p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])
                                            
                                            if p0 == nothing
                                                int = false
                                            end
                                            #--------------
                                        end

                                    else
                                        bo[i].ver1[ind_v2,1] = xtmp
                                        bo[i].ver1[ind_v2,2] = ztmp

                                        #--------------
                                        p2= [bo[i].ver1[ind_v2,1], bo[i].ver1[ind_v2,2]]
                                        p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])
                                        
                                        if p0 == nothing
                                            int = false
                                        end
                                        #--------------
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if !int
                #idpolyorg = deepcopy(idpoly)
                for k in collect(1:length(borg))
                    borg[k].ver1[:] = deepcopy(bo[k].ver1[:])                       
                end
            end
        end
    end
    return nothing #idpolyorg        #idpoly: array constituted by other arrays (length = n° polygons) containing the index
                            #of the first vertex of all sides crossing another side
end

#####################################################################################################

"""
$(TYPEDSIGNATURES)

Function to check if there are intersections beetween sides of a polygon in respect to neighbour polygons, returning a Boolean value (i.e. true or false).
"""
function intersectpairpoly(bo::Vector{BodySegments2D})
    
    nbo = shuffle(collect(1:length(bo)))

    # No check if only one polygon
    if length(nbo)==1
        return false
    end    

    numpoly = copy(nbo)	
    
    for i in nbo

        deleteat!(numpoly, findall(x->x==i,numpoly))
                
    	if isempty(numpoly)
            break
	end

        shuffle!(numpoly)
        nedge1 = length(bo[i].ver1[:,1])
        polyarray1 = Array{Array{Float64,1},1}(undef,(nedge1+1))

        for n=1:nedge1
            polyarray1[n] = bo[i].ver1[n,:]
            if n==nedge1
                polyarray1[end] = bo[i].ver1[1,:]
            end
        end
        
        # define polygon1 for LibGEOS
        poly1 = LibGEOS.Polygon([polyarray1])
        
        for j in numpoly
            
            nedge2 = length(bo[j].ver1[:,1])
            polyarray2 = Array{Array{Float64,1},1}(undef,(nedge2+1))
            
            for n=1:nedge2
                polyarray2[n] = bo[j].ver1[n,:]
                if n==nedge2
                    polyarray2[end] = bo[j].ver1[1,:]
                end
            end
            # define polygon2 for LibGEOS
            poly2 = LibGEOS.Polygon([polyarray2])
            
            # assess polygons overlapping
            if LibGEOS.overlaps(poly1,poly2)
                return true
            end
        end
    end
    return false
end


#####################################################################################################

"""
$(TYPEDSIGNATURES)

Function to check if there are intersections beetween sides of a polygon in respect to neighbour polygons, trying to shift polygon vertices in order to avoid intersection. This algorithm shifts the polygon vertices in agreement with the mathematical properties of the Hamiltonian Dynamics.
"""
function verpolyallshift!(borg::Vector{BodySegments2D},
                          blow::Union{Vector{BodySegments2D},Nothing}=nothing,
                          bup::Union{Vector{BodySegments2D},Nothing}=nothing)

    @assert typeof(blow)==typeof(bup)

    # No check if only one polygon
    if length(borg)==1
        return nothing
    end
    
    # idpolyorg = [Int64[] for a in 1:length(borg)]
    fix = true
    countfix = 0
    sign1 = rand(2:5)*length(borg)
    sign2 = rand(2:5)*length(borg)
    
    while fix

        if countfix == sign1
            break
        end

        countfix += 1
        # idpoly = deepcopy(idpolyorg)
        fix = false
        bo = deepcopy(borg)
        nbo = shuffle(collect(1:length(bo)))
        int = true
        count = 0
        
        while int

            if count == sign2
                fix = true
                break
            end
            
            count += 1 
            int = false
            numpoly = copy(nbo)
            
            for i in nbo
                num_vert1 = length(bo[i].ver1[:,1])
                deleteat!(numpoly, findall(x->x==i,numpoly))
                
    	        if isempty(numpoly)
	            break
	        end
                
                shuffle!(numpoly)
                
                for ind_v1=1:num_vert1
                    
                    if ind_v1 < num_vert1
                        p1= [bo[i].ver1[ind_v1,1], bo[i].ver1[ind_v1,2]]
                        q1= [bo[i].ver1[ind_v1+1,1], bo[i].ver1[ind_v1+1,2]]
                    else
                        p1= [bo[i].ver1[ind_v1,1], bo[i].ver1[ind_v1,2]]
                        q1= [bo[i].ver1[1,1], bo[i].ver1[1,2]]
                    end
                    
                    for j in numpoly
                        
                        num_vert2 = length(bo[j].ver1[:,1])
                        
                        for ind_v2=1:num_vert2 #taking into account last segment vertice -> first vertice
                            
                            p2= [bo[j].ver1[ind_v2,1], bo[j].ver1[ind_v2,2]]
                            
                            if ind_v2==num_vert2
                                
                                q2= [bo[j].ver1[1,1], bo[j].ver1[1,2]]
                                
                                doInternalp2 = isInternal(bo[i],p2[1],p2[end])
                                doInternalq2 = isInternal(bo[i],q2[1],q2[end])
                                
                                p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])

                                if p0 != nothing
                                    int = true
                                end
                                
                                if (doInternalp2 || doInternalq2) && p0 != nothing  
                                    
                                    if doInternalp2
                                        
                                        dist1 = sqrt(((p2[1]-q2[1])^2)+((p2[end]-q2[end])^2))
                                        dist2 = sqrt(((p0[1]-p2[1])^2)+((p0[end]-p2[end])^2))   
                                        ratio = dist2/dist1
                                        Δx = p0[1]-q2[1]
                                        Δz = p0[end]-q2[end]
                                        xtmp = p0[1] - (Δx*ratio)
                                        ztmp = p0[end] - (Δz*ratio)
                                        
                                        if blow != nothing && bup != nothing
                                            
                                            if (blow[j].ver1[ind_v2,1] <= xtmp <= bup[j].ver1[ind_v2,1])&&
                                                (blow[j].ver1[ind_v2,2] <= ztmp <= bup[j].ver1[ind_v2,2])
                                                
                                                bo[j].ver1[ind_v2,1] = xtmp
                                                bo[j].ver1[ind_v2,2] = ztmp

                                                #--------------
                                                p2= [bo[j].ver1[ind_v2,1], bo[j].ver1[ind_v2,2]]
                                                p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])

                                                if p0 == nothing
                                                    int = false
                                                end
                                                #--------------
                                            end
                                            
                                        else
                                            bo[j].ver1[ind_v2,1] = xtmp
                                            bo[j].ver1[ind_v2,2] = ztmp
                                            
                                            #--------------
                                            p2= [bo[j].ver1[ind_v2,1], bo[j].ver1[ind_v2,2]]
                                            p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])

                                            if p0 == nothing
                                                int = false
                                            end
                                            #--------------
                                        end
                                        
                                    elseif doInternalq2
                                        
                                        dist1 = sqrt(((p2[1]-q2[1])^2)+((p2[end]-q2[end])^2))
                                        dist2 = sqrt(((p0[1]-q2[1])^2)+((p0[end]-q2[end])^2))   
                                        ratio = dist2/dist1
                                        Δx = p0[1]-p2[1]
                                        Δz = p0[end]-p2[end]
                                        xtmp = p0[1] - (Δx*ratio)
                                        ztmp = p0[end] - (Δz*ratio)
                                            
                                        if blow != nothing && bup != nothing
                                            
                                            if (blow[j].ver1[1,1] <= xtmp <= bup[j].ver1[1,1])&&
                                                (blow[j].ver1[1,2] <= ztmp <= bup[j].ver1[1,2])

                                                bo[j].ver1[1,1] = xtmp
                                                bo[j].ver1[1,2] = ztmp

                                                #--------------
                                                q2= [bo[j].ver1[1,1], bo[j].ver1[1,2]]
                                                p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])

                                                if p0 == nothing
                                                    int = false
                                                end
                                                #--------------
                                            end
                                            
                                        else
                                            bo[j].ver1[1,1] = xtmp
                                            bo[j].ver1[1,2] = ztmp

                                            #--------------
                                            q2= [bo[j].ver1[1,1], bo[j].ver1[1,2]]
                                            p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])

                                            if p0 == nothing
                                                int = false
                                            end
                                            #--------------
                                        end
                                    end                              
                                end
                            else
                                
                                q2= [bo[j].ver1[ind_v2+1,1], bo[j].ver1[ind_v2+1,2]]
                                
                                doInternalp2 = isInternal(bo[i],p2[1],p2[end])
                                doInternalq2 = isInternal(bo[i],q2[1],q2[end])
                                
                                p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])

                                if p0 != nothing
                                    int = true
                                end

                                if (doInternalp2 || doInternalq2) && p0 != nothing  
                                                                        
                                    if doInternalp2
                                        
                                        dist1 = sqrt(((p2[1]-q2[1])^2)+((p2[end]-q2[end])^2))
                                        dist2 = sqrt(((p0[1]-p2[1])^2)+((p0[end]-p2[end])^2))   
                                        ratio = dist2/dist1
                                        Δx = p0[1]-q2[1]
                                        Δz = p0[end]-q2[end]
                                        xtmp = p0[1] - (Δx*ratio)
                                        ztmp = p0[end] - (Δz*ratio)
                                            
                                        if blow != nothing && bup != nothing
                                            
                                            if (blow[j].ver1[ind_v2,1] <= xtmp <= bup[j].ver1[ind_v2,1])&&
                                                (blow[j].ver1[ind_v2,2] <= ztmp <= bup[j].ver1[ind_v2,2])
                                                
                                                bo[j].ver1[ind_v2,1] = xtmp
                                                bo[j].ver1[ind_v2,2] = ztmp

                                                #--------------
                                                p2= [bo[j].ver1[ind_v2,1], bo[j].ver1[ind_v2,2]]
                                                p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])

                                                if p0 == nothing
                                                    int = false
                                                end
                                                #--------------
                                            end
                                            
                                        else
                                            bo[j].ver1[ind_v2,1] = xtmp
                                            bo[j].ver1[ind_v2,2] = ztmp

                                            #--------------
                                            p2= [bo[j].ver1[ind_v2,1], bo[j].ver1[ind_v2,2]]
                                            p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])

                                            if p0 == nothing
                                                int = false
                                            end
                                            #--------------
                                        end
                                        
                                    elseif doInternalq2
                                        
                                        dist1 = sqrt(((p2[1]-q2[1])^2)+((p2[end]-q2[end])^2))
                                        dist2 = sqrt(((p0[1]-q2[1])^2)+((p0[end]-q2[end])^2))   
                                        ratio = dist2/dist1
                                        Δx = p0[1]-p2[1]
                                        Δz = p0[end]-p2[end]
                                        xtmp = p0[1] - (Δx*ratio)
                                        ztmp = p0[end] - (Δz*ratio)
                                            
                                        if blow != nothing && bup != nothing
                                            
                                            if (blow[j].ver1[ind_v2+1,1] <= xtmp <= bup[j].ver1[ind_v2+1,1])&&
                                                (blow[j].ver1[ind_v2+1,2] <= ztmp <= bup[j].ver1[ind_v2+1,2])
                                                
                                                bo[j].ver1[ind_v2+1,1] = xtmp
                                                bo[j].ver1[ind_v2+1,2] = ztmp

                                                #--------------
                                                q2= [bo[j].ver1[ind_v2+1,1], bo[j].ver1[ind_v2+1,2]]
                                                p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])

                                                if p0 == nothing
                                                    int = false
                                                end
                                                #--------------
                                            end
                                            
                                        else
                                            bo[j].ver1[ind_v2+1,1] = xtmp
                                            bo[j].ver1[ind_v2+1,2] = ztmp

                                            #--------------
                                            q2= [bo[j].ver1[ind_v2+1,1], bo[j].ver1[ind_v2+1,2]]
                                            p0 = Inter2Segm(p1[1],p1[end],q1[1],q1[end],p2[1],p2[end],q2[1],q2[end])

                                            if p0 == nothing
                                                int = false
                                            end
                                            #--------------
                                        end
                                    end                                                                                             
                                end
                            end
                        end
                    end
                end
            end
            if !int
                # idpolyorg = deepcopy(idpoly)
                for k in collect(1:length(borg))
                    borg[k].ver1[:] = deepcopy(bo[k].ver1[:])                       
                end
            end
        end
    end
    
    return nothing #idpolyorg           #idpoly: array containing other arrays (length = n° polygons) containing the index                                 
                               #of the first vertex of all sides crossing another side
end

###################################################################################

