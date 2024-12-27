###################################################

"""
    Calculate area of a set of polygons given a `PolygBodies2D` structure.
    Returns a vector where each element is the area of one of the polygons.
"""
function calcareamanypoly(pbody::PolygBodies2D) 
    nbo = length(pbody.bo)
    area = Vector{Float64}(undef,nbo)
    for b=1:nbo
        # body segments
        bs = pbody.bo[b]
        # area
        area[b] = calcareapoly(bs)
    end
    return area
end

###################################################

"""
    Calculate the area of a single polygon given a `BodySegments2D` structure.
    Returns the value of the area.
"""
function calcareapoly(bs::BodySegments2D) 
    # https://mathworld.wolfram.com/PolygonArea.html
    # https://stackoverflow.com/a/10298685
    encarea2=0.0 # twice the area
    for ise=1:bs.nsegm
        x1 = bs.ver1[ise,1]
        z1 = bs.ver1[ise,2] 
        x2 = bs.ver2[ise,1]
        z2 = bs.ver2[ise,2]
        encarea2 += (x2-x1)*(z2+z1)
    end
    # final area of the polygon
    # abs because area is
    #      positive if polygon is counterclockwise
    #      negative if polygon is sclockwise
    area = abs(encarea2 / 2.0) 
    return area
end

###################################################
