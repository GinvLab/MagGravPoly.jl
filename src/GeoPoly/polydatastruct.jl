


##########################################
"""
$(TYPEDEF)

Structure containing the segments of a polygonal body.
To create an instance a set of indices have to be passed on.

# Fields 

$(TYPEDFIELDS)

"""
struct BodySegments2D
    "(x,y) for first set of vertices (beginning of segments)"
    ver1::SubArray 
    "(x,y) for second set of vertices (end of segments)"
    ver2::SubArray
    "total number of segments"
    nsegm::Integer

    function BodySegments2D(idx1::Vector{<:Integer},vertices::AbstractArray{<:Real,2})
        
        # @assert ndims(ver1)==2
        @assert size(vertices,2)==2
        # circular shift to get second set of indices
        idl = length(idx1)
        ids = vcat([2:idl;],1)
        idx2 = view(idx1,ids)
        ## first set of vertices
        ver1 = view(vertices,idx1,:)
        ## second set of vertices
        ver2 = view(vertices,idx2,:)
        nsegm = size(ver1,1)
        return new(ver1,ver2,nsegm)
    end    
end

##########################################
"""
$(TYPEDEF)

Structure containing a set of polygonal bodies (described by their segments and all vertices).
To create an instance, input an array of vectors of indices (of vertices) for each body and the array of all the vertices.

# Fields 

$(TYPEDFIELDS)

"""
struct PolygBodies2D
    "array of bodies defined by their vertices"
    bo::Vector{BodySegments2D}
    "number of polygonal bodies"
    nbo::Integer
    "array of all vertices for all bodies"
    allvert::AbstractArray{<:Real,2}
    "array indices relating to the vertices in allvert"
    bodyindices::Vector{<:Vector{<:Integer}}

    function PolygBodies2D(bodyindices::Vector{<:Vector{<:Integer}},allvert::AbstractArray{<:Real,2})

        @assert all((length.(bodyindices)).>=3)
        N=length(bodyindices)
        #Check on bodyindices
        # if N > 1
        #     indbo = let
        #         num = length(bodyindices[1])
        #         for i=2:N
        #             N2 = collect(1:N)
        #             deleteat!(N2, findall(x->x>=i,N2))
        #             for j in N2
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

        ## Get all the unique body indices
        #allind = vcat(bodyindices...)
        allind = collect(Iterators.flatten(bodyindices))
        unidxs = unique(allind)

        @assert length(unidxs)==length(allvert[:,1])
        @assert size(allvert,2)==2
        @assert size(allvert,1)>=3 # it has to be at least a triangle...
        # make a copy to make sure that pointers point only to internal stuff
        allvert_copy = copy(allvert)
        bodyindices_copy = deepcopy(bodyindices)
        #
        bo = Vector{BodySegments2D}(undef,N)
        for i=1:N
            bo[i] = BodySegments2D(bodyindices_copy[i],allvert_copy)
        end
        # # induced magnetization
        # Jind = MagnetizVector(mod=modind,Ideg=Idegind,Ddeg=Ddegind)
        # # remnant magnetization
        # Jrem = MagnetizVector(mod=modrem,Ideg=Idegrem,Ddeg=Ddegrem)
        return new(bo,N,allvert_copy,bodyindices_copy)   
    end    
end


###############################################

"""
$(TYPEDSIGNATURES)

Julia structure defining a topography as set of segments characterized by x and z coordinates and both angular coefficient and intercept.  

# Fields 

$(TYPEDFIELDS)

"""
struct TopoEdges
    
    "all x couples of topography coordinates"
    verx::Array{Array{Float64,1},1}
    "all z couples of topography coordinates"
    verz::Array{Array{Float64,1},1}
    "angular coefficients and segment intersepts of each coordinates couple"
    mq::Array{Array{Float64,1},1}
    "2D array of (x,z) coordinates for topography"
    xz::Array{Float64,2}

    function TopoEdges(topo::Array{Float64,2})
        # @assert ndims(topo)==2
        @assert size(topo,2)==2
        # definition array of ver and mq
        verx = Array{Array{Float64,1}}(undef,length(topo[:,1])-1)
        verz = Array{Array{Float64,1}}(undef,length(topo[:,1])-1)
        mq = Array{Array{Float64,1}}(undef,length(topo[:,1])-1)

        for i=1:length(topo[:,1])-1
            verx[i] = [topo[i,1],topo[i+1,1]]
            verz[i] = [topo[i,2],topo[i+1,2]]
            Δx = topo[i+1,1]-topo[i,1]
            Δz = topo[i+1,2]-topo[i,2]
            m = Δz/Δx
            q = topo[i,2] - m*topo[i,1]
            mq[i] = [m,q]
        end

        xz = copy(topo)
        return new(verx,verz,mq,xz)
    end    
end

##########################################
