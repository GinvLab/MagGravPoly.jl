#######################################################################################
"""
$(TYPEDSIGNATURES)

Check whether the polygonal body has segments ordered anticlockwise.
"""
function checkanticlockwiseorder(body::BodySegments2D)
    ## https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
    ## https://www.element84.com/blog/determining-the-winding-of-a-polygon-given-as-a-set-of-ordered-points
    #
    # Check direction (anti)clockwise for a reference
    #   system like the following:
    #
    #   z
    #  /\ 
    #  |      2
    #  |   1     3
    #  |      4
    #  |
    #  -------------> x
    #
    encarea2=0.0 # twice the area
    for ise=1:body.nsegm
        x1 = body.ver1[ise,1]
        z1 = body.ver1[ise,2] 
        x2 = body.ver2[ise,1]
        z2 = body.ver2[ise,2]
        encarea2 += (x2-x1)*(z2+z1)
    end
    #
    # anticlockwise -> encarea2 < 0.0
    # clockwise -> encarea2 > 0.0
    if encarea2<0.0
        anticlockw=true
    else
        anticlockw=false
    end
    #
    # The reference system for the magnetic anomaly functions
    #   is reversed in z:
    #
    #  -------------> x
    #  |
    #  |      4
    #  |   1     3
    #  |      2
    #  \/
    #  z
    #
    # so, consequently, we flip the direction of
    # clockwise/anticlockwise:   !(anticlockw)   
    return !anticlockw    
end


#######################################################################################
#=
"""
$(TYPEDSIGNATURES)

Function to check all polygons' vertices and bodyindices in order to make sure that vertices in 
 common among polygons appear only once.
"""
function checkbodyindices(bodyindicesorg::Vector{<:Vector{<:Integer}},verticesorg::AbstractArray{<:Real,2})

    #-------------------------------------
    function checktype(arr::Vector{Int64})
        type = "Up"
        for i=1:length(arr)-1
            if arr[i] >= arr[i+1]
                type = "Down"
                return type
            end
        end
        return type
    end
    #-------------------------------------
    
    vertices = deepcopy(verticesorg)
    bodyindices = deepcopy(bodyindicesorg)
    nbo = length(bodyindices)

    ## All checks for more than one polygon
    for i=1:nbo
        type=checktype(bodyindices[i])

        if type == "Down"
            error("Bodyindices must contain values in increasing order. Perhaps check on bodyindices already performed. Aborting!")
        end

        if i<nbo
            if (intersect(bodyindices[i],bodyindices[i+1])) != []
                error("Bodyindices for each polygon must contain different elements. Perhaps check on bodyindices already performed. Aborting!")
            end
        end

        # No other checks if only one polygon
        if nbo == 1
            pbody = PolygBodies2D(bodyindices,vertices)
            check=checkanticlockwiseorder(pbody.bo[i])
            if !check
                bodyindices[i] = reverse(bodyindices[i],dims = 1)
            end
            return bodyindices,vertices
        end
    end

    idarray=collect(1:length(vertices[:,1]))
    minidx = Int64[]
    idx = Array{Int64}(undef, 0, 2)
    idallvertfix = Array{Int64,2}(undef,0,2)
    i=1

    while i != length(idarray)

        sign = [vertices[i,:].==vertices[j,:] for j in collect(1:length(vertices[:,1]))]
        id = findall(x->x==[1,1],sign)

        if length(id)>=2

            minid = minimum(id)
            idallvert = copy(id)
            filter!(x->x!=minid,idallvert)

            for a in id

                for b in collect(1:nbo)

                    if minidx == []
                        minidbody=findall(x->x==minid,bodyindices[b])

                        if minidbody != []
                            append!(minidx,[b,minidbody[1]])
                        end
                    end

                    idbody=findall(x->x==a,bodyindices[b])
                    idbodyorg = copy(idbody)

                    if b>1 && idbody != []
                        for l=1:(b-1)
                            idbodyorg[1] += length(bodyindices[l])
                        end
                    end

                    if (idbody != []) && (minidx != []) && (idbodyorg[1] != minidx[2])
                        idx = vcat(idx,[b idbody[1]])
                        idallvertfix = vcat(idallvertfix,[b idbody[1]])
                    end
                end
            end

            for c in collect(1:length(idx[:,1]))

                bodyindices[idx[c,1]][idx[c,2]]=copy(bodyindices[minidx[1]][minidx[2]])
                idallvertfixc=copy(idallvertfix)

                if (ind = findall(x->x!=idx[c,1],idallvertfixc[:,1])) != []
                    idallvertfixc=idallvertfixc[setdiff(1:end, ind), :]
                end

                if idx[c,2] != length(bodyindices[idx[c,1]])

                    for d in collect(idx[c,2]+1:length(bodyindices[idx[c,1]]))

                        if findall(x->x==1,d .==idallvertfixc[:,2]) == []
                            bodyindices[idx[c,1]][d] -= 1
                        end
                    end
                end
            end

            if maximum(idx[:,1])!=nbo
                for k=1:nbo-maximum(idx[:,1])
                    bodyindices[maximum(idx[:,1])+k][1:end] .-= 1
                end
            end

            minidx = Int64[]
            idx = Array{Int64}(undef, 0, 2)
            vertices = vertices[setdiff(1:end, idallvert), :]
            idarray = collect(1:length(vertices[:,1]))

        end
        i += 1
    end

    # Checking anti-clockwise orientation of polygons
    pbody = PolygBodies2D(bodyindices,vertices)
    
    for l=1:nbo 
        check=checkanticlockwiseorder(pbody.bo[l])
        if !check
            bodyindices[l] = reverse(bodyindices[l],dims = 1)
        end
    end

    return bodyindices,vertices
end
=#
##############

"""
$(TYPEDSIGNATURES)

Function to check for consistency between polygon vertex indices (`bodyindices`) and vertex coordinates (`vertices`). It performs the following checks and corrections:

1. **Monotonic ordering** – Each polygon's vertex indices must be strictly increasing.
2. **Index uniqueness** – No two polygons should share the same vertex indices.
3. **Duplicate vertex merging** – Identical vertices are merged, and indices are updated.
4. **Counterclockwise orientation** – All polygons are reoriented to counterclockwise order.

# Returns
A tuple `(bodyindices, vertices)` where:
- `bodyindices` are the corrected vertex index lists for each polygon.
- `vertices` is the updated vertex matrix with duplicates removed.
"""
function checkbodyindices(bodyindices::Vector{<:Vector{<:Integer}},
                          vertices::AbstractMatrix{<:Real})

    # -------------------------------
    # Internal helper: check if array is strictly increasing
    # -------------------------------
    check_increasing(arr::Vector{Int}) = all(diff(arr) .> 0)

    # Make deep copies to avoid mutating the original input
    bodyindices = deepcopy(bodyindices)
    vertices    = deepcopy(vertices)
    nbo         = length(bodyindices)

    # -------------------------------
    # Step 1: Check each polygon’s index ordering and uniqueness
    # -------------------------------
    for (i, inds) in enumerate(bodyindices)
        # Ensure indices are sorted in ascending order
        if !check_increasing(inds)
            error("Indices in bodyindices[$i] must be strictly increasing. Aborting!")
        end

        # Ensure that consecutive polygons do not share vertex indices
        if i < nbo && !isempty(intersect(bodyindices[i], bodyindices[i + 1]))
            error("Polygons $i and $(i + 1) share vertex indices, which is not allowed. Aborting!")
        end
    end

    # -------------------------------
    # Step 2: Handle the simple case (single polygon)
    # -------------------------------
    if nbo == 1
        pbody = PolygBodies2D(bodyindices, vertices)
        if !checkanticlockwiseorder(pbody.bo[1])
            bodyindices[1] = reverse(bodyindices[1])
        end
        return bodyindices, vertices
    end

    # -------------------------------
    # Step 3: Merge duplicate vertices
    # -------------------------------
    i = 1
    while i <= size(vertices, 1)
        # Find all vertices identical to vertex i
        duplicates = findall(j -> all(vertices[i, :] .== vertices[j, :]),
                             1:size(vertices, 1))

        if length(duplicates) > 1
            # Keep the first as reference
            ref = first(duplicates)
            dupes = setdiff(duplicates, [ref])

            # Update all polygons to replace duplicates with the reference index
            for b in eachindex(bodyindices)
                bodyindices[b] = [x in dupes ? ref : x for x in bodyindices[b]]
            end

            # Remove duplicate vertices from the matrix
            vertices = vertices[setdiff(1:end, dupes), :]

            # Adjust indices after removal of vertices
            for b in eachindex(bodyindices)
                bodyindices[b] = [x - count(d -> d < x, dupes) for x in bodyindices[b]]
            end
        end
        i += 1
    end

    # -------------------------------
    # Step 4: Ensure polygons are counterclockwise
    # -------------------------------
    pbody = PolygBodies2D(bodyindices, vertices)
    for (i, bo) in enumerate(pbody.bo)
        if !checkanticlockwiseorder(bo)
            bodyindices[i] = reverse(bodyindices[i])
        end
    end

    return bodyindices, vertices
end

#######################################################

