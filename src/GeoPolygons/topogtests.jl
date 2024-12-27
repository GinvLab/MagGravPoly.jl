#########################

"""
$(TYPEDSIGNATURES)

Function to check if there is some polygon crossing the topography, returning a Boolean value (i.e. true or false).  
"""
function checktopo(topo::TopoEdges,bo::Vector{BodySegments2D}) 

    ibo = shuffle(collect(1:length(bo)))
    nbo = length(bo)

    ##======================================================
    ## construct a polygon from the topography array (which
    ##    is only an open polyline)
    nedgetopo = size(topo.xz,1)+3 
    topoarray = Array{Array{Float64,1},1}(undef,nedgetopo)   

    # highest among vertex and topography
    zmin_topo = minimum(topo.xz[:,2]) ## z is pointing downward!
    zmin_vert = minimum([minimum(bo[i].ver1[:,2]) for i=1:nbo]) ## z is pointing downward!
    zminsc = 1.5 * min(zmin_topo,zmin_vert) ## z is pointing downward!
 
    for n=1:nedgetopo
        if n==nedgetopo
            topoarray[end] = topo.xz[1,:]
        elseif n==nedgetopo-1
            topoarray[n] = [topo.xz[1,1],zminsc]
        elseif n==nedgetopo-2
            topoarray[n] = [topo.xz[n-1,1],zminsc]
        else
            topoarray[n] = topo.xz[n,:]
        end
    end
    
    polytopo = LibGEOS.Polygon([topoarray])

    ##======================================================
    ## create polygons usable by LibGEOS from the
    ##   polygonal bodies
    for i in ibo
        
    	nedge = length(bo[i].ver1[:,1])
        polyarray = Array{Array{Float64,1},1}(undef,(nedge+1))   
        
        for n=1:nedge
            polyarray[n] = bo[i].ver1[n,:]
            if n==nedge
                polyarray[end] = bo[i].ver1[1,:]
            end
        end
        
        poly = LibGEOS.Polygon([polyarray])

        if any(any.([isnan.(to) for to in topoarray])) || any(any.([isnan.(po) for po in polyarray]))
            println("\n### checktopo(): any(any.([isnan.(to) for to in topoarray])) || any(any.([isnan.(po) for po in polyarray]))  ###\n")
        end

        if LibGEOS.overlaps(polytopo,poly)
            return true
        end
    end
    return false
end
    
#####################################

"""
$(TYPEDSIGNATURES)

Function to check if there is some polygon crossing the topography, shifting the polygon vertices in order to avoid intersection. This algorithm shifts the polygon vertices in agreement with the mathematical properties of the Hamiltonian Dynamics.  
"""
function vertoposhift!(topo::TopoEdges,bo::Vector{BodySegments2D},
                       blow::Union{Vector{BodySegments2D},Nothing}=nothing,
                       bup::Union{Vector{BodySegments2D},Nothing}=nothing)

    @assert typeof(blow)==typeof(bup)
    nbo = length(bo)
    #idtopofix = [Int64[] for a in 1:nbo]

    for i=1:nbo

        for j=1:length(topo.verx)

            for k=1:length(bo[i].ver1[:,1])
                
                ind = findall(x -> topo.verx[j][1]<= x <=topo.verx[j][2], bo[i].ver1[k,1])
                
                if ind!=[]
                    
                    ztopo = topo.mq[j][1]*bo[i].ver1[k,1]+topo.mq[j][2]
                    
                    if bo[i].ver1[k,2] < ztopo

                        if blow != nothing && bup != nothing
                            
                            if (blow[i].ver1[k,2]<=ztopo<=bup[i].ver1[k,2])
                                
                                bo[i].ver1[k,2] = ztopo #2*ztopo - bo[i].ver1[k,2]
                            end
                            
                        else
                            bo[i].ver1[k,2] = ztopo #2*ztopo - bo[i].ver1[k,2]       
                        end     
                        #append!(idtopofix[i],k)
                    end
                end
            end
        end 
    end
    return nothing #idtopofix
end

######################################################

"""
$(TYPEDSIGNATURES)

Function to check if the topography extends both righward and leftward beyond the outermost polygons by a user-defined percentage (> 25%).  
"""
function checkmodelizdim(topo::TopoEdges,bo::Vector{BodySegments2D},perc::Real)

    if perc < 25.0
        error("The minimum value of the 'perc' input is 25.0. Aborting")
    end
    
    maxallx = -Inf
    minallx = Inf

    for i=1:length(bo)
        maxx = maximum(bo[i].ver1[:,1])
        minx = minimum(bo[i].ver1[:,1])
        if maxx > maxallx
            maxallx = copy(maxx)
        end
        if minx < minallx
            minallx = copy(minx)
        end
    end

    topomaxx = maximum(topo.xz[:,1])
    topominx = minimum(topo.xz[:,1])
    Δx2topo = abs(topomaxx-topominx)/2
    
    if maxallx < topomaxx && minallx > topominx

        Δx2maxpoly = abs(topomaxx-maxallx)
        Δx2minpoly = abs(topominx-minallx)  
        percmaxx = (Δx2maxpoly/Δx2topo)*100.0
        percminx = (Δx2minpoly/Δx2topo)*100.0
        
        if percmaxx < perc && percminx < perc
            ext2maxx = floor((((perc*Δx2topo)/100.0) - Δx2maxpoly),digits=2)
            ext2minx = floor((((perc*Δx2topo)/100.0) - Δx2minpoly),digits=2)
            error("Topography must be extended at least $ext2maxx units to the right and $ext2minx units to the left. Aborting")
        elseif percmaxx < perc && percminx >= perc
            ext2maxx = floor((((perc*Δx2topo)/100.0) - Δx2maxpoly),digits=2)
            error("Topography must be extended at least $ext2maxx units to the right. Aborting")
        elseif percmaxx >= perc && percminx < perc
            ext2minx = floor((((perc*Δx2topo)/100.0) - Δx2minpoly),digits=2)
            error("Topography must be extended at least $ext2minx units to the left. Aborting")
        end

    elseif maxallx < topomaxx && minallx <= topominx

        Δx2maxpoly = abs(topomaxx-maxallx)
        percmaxx = (Δx2maxpoly/Δx2topo)*100.0   
        topo2minx = abs(topominx-minallx)
        ext2minx = floor((((perc*Δx2topo)/100.0) + topo2minx),digits=2)

        if percmaxx < perc
            ext2maxx = floor((((perc*Δx2topo)/100.0) - Δx2maxpoly),digits=2)
            error("Topography must be extended at least $ext2maxx units to the right and $ext2minx units to the left. Aborting")
        else
            error("Topography must be extended at least $ext2minx units to the left. Aborting")
        end
        
    elseif maxallx >= topomaxx && minallx > topominx

        Δx2minpoly = abs(topominx-minallx)
        percminx = (Δx2minpoly/Δx2topo)*100.0
        topo2maxx = abs(topomaxx-maxallx)
        ext2maxx = floor((((perc*Δx2topo)/100.0) + topo2maxx),digits=2)

        if percminx < perc
            ext2minx = floor((((perc*Δx2topo)/100.0) - Δx2minpoly),digits=2)
            error("Topography must be extended at least $ext2maxx units to the right and $ext2minx units to the left. Aborting")
        else
            error("Topography must be extended at least $ext2maxx units to the right. Aborting")
        end
        
    elseif maxallx >= topomaxx && minallx <= topominx

        topo2minx = abs(topominx-minallx)
        topo2maxx = abs(topomaxx-maxallx)
        ext2minx = floor((((perc*Δx2topo)/100.0) + topo2minx),digits=2)
        ext2maxx = floor((((perc*Δx2topo)/100.0) + topo2maxx),digits=2)
        error("Topography must be extended at least $ext2maxx units to the right and $ext2minx units to the left. Aborting")
        
    end
    
    return nothing
end


##################################################
############################################################

