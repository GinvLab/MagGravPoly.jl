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

##############
