#--------------------------------------------------------#
#---------------------PACKAGE-TESTS----------------------#
#--------------------------------------------------------#


using MagGravPoly.GeoPolygons

########################################

## Single polygon test
function poly_assembly1()

    vertices  = [350.0 500.0;
                 650.0 500.0;
                 800.0 350.0;
                 650.0 200.0;
	         350.0 200.0;
	         200.0 350.0]

    ind = collect(1:6)
    bodyindices = [ind]
    pbody = PolygBodies2D(bodyindices,vertices)

    ##----------------------------------------------------------
    # Polygon reverse
    reverse!(pbody.bodyindices[1])
    
    ## list of ver1
    ver1 = [200.0 350.0;
            350.0 200.0;
            650.0 200.0;
            800.0 350.0;
            650.0 500.0;
            350.0 500.0]
    
    ## list of ver2
    ver2 = [350.0 200.0;
            650.0 200.0;
            800.0 350.0;
            650.0 500.0;
            350.0 500.0;
            200.0 350.0]
    
    ##----------------------------------------------
    # Check if lists of starting and ending segment vertices mach (test passed)
    same1 = all(ver1.==pbody.bo[1].ver1)
    same2 = all(ver2.==pbody.bo[1].ver2)
    
    if same1 && same2
        return true
    elseif same1 && !same2
        printstyled("\nField `ver2` in `BodySegments2D` containing wrong list of vertices. Test failed! \n", bold=true,color=:red)
        return false
    elseif !same1 && same2
        printstyled("\nField `ver1` in `BodySegments2D` containing wrong list of vertices. Test failed! \n", bold=true,color=:red)
        return false
    end
end


#######################################################################

## Two polygons test
function poly_assembly2()

    vertices  = [350.0 500.0;
                 650.0 500.0;
                 800.0 350.0;
                 650.0 200.0;
	         350.0 200.0;
	         200.0 350.0]
    
    ind1 = collect(1:3)
    ind2 = collect(4:6)
    bodyindices = [ind1,ind2]
    pbody = PolygBodies2D(bodyindices,vertices)
    
    ##----------------------------------------------------------
    # Polygon reverse
    id = rand(1:3)
    
    # Multiple cases
    if id == 1
        reverse!(pbody.bodyindices[id])

        ## lists of ver1 (poly 1 & poly 2)
        ver11 = [800.0 350.0;
                 650.0 500.0;
                 350.0 500.0]

        ver12 = [650.0 200.0;
                 350.0 200.0;
                 200.0 350.0]
        
        ## lists of ver2 (poly 1 & poly 2)
        ver21 = [650.0 500.0;
                 350.0 500.0;
                 800.0 350.0]

        ver22 = [350.0 200.0;
                 200.0 350.0;
                 650.0 200.0]
        
    elseif id == 2
        reverse!(pbody.bodyindices[id])

        ## lists of ver1 (poly 1 & poly 2)
        ver11 = [350.0 500.0;
                 650.0 500.0;
                 800.0 350.0]

        ver12 = [200.0 350.0;
                 350.0 200.0;
                 650.0 200.0]
        
        ## lists of ver2 (poly 1 & poly 2)
        ver21 = [650.0 500.0;
                 800.0 350.0;
                 350.0 500.0]

        ver22 = [350.0 200.0;
                 650.0 200.0;
                 200.0 350.0]
        
    elseif id == 3
        reverse!(pbody.bodyindices[1])
        reverse!(pbody.bodyindices[2])

        ## lists of ver1 (poly 1 & poly 2)
        ver11 = [800.0 350.0;
                 650.0 500.0;
                 350.0 500.0]

        ver12 = [200.0 350.0;
                 350.0 200.0;
                 650.0 200.0]

        ## lists of ver2 (poly 1 & poly 2)
        ver21 = [650.0 500.0;
                 350.0 500.0;
                 800.0 350.0]

        ver22 = [350.0 200.0;
                 650.0 200.0;
                 200.0 350.0]
        
    end
    
    
    ##----------------------------------------------
    # Check if lists of starting and ending segment vertices mach (test passed)
    same11 = all(ver11.==pbody.bo[1].ver1)
    same12 = all(ver12.==pbody.bo[2].ver1)
    same21 = all(ver21.==pbody.bo[1].ver2)
    same22 = all(ver22.==pbody.bo[2].ver2)
    
    if same11 && same12 && same21 && same22
        return true
    else
        if !same11
            printstyled("\nField `ver1` in `BodySegments2D` of 1st polygon containing wrong list of vertices. Test failed! \n", bold=true,color=:red)
        elseif !same12
            printstyled("\nField `ver1` in `BodySegments2D` of 2nd polygon containing wrong list of vertices. Test failed! \n", bold=true,color=:red)
        elseif !same21
            printstyled("\nField `ver2` in `BodySegments2D` of 1st polygon containing wrong list of vertices. Test failed! \n", bold=true,color=:red)
        elseif !same22
            printstyled("\nField `ver2` in `BodySegments2D` of 2nd polygon containing wrong list of vertices. Test failed! \n", bold=true,color=:red)    
        end
        return false
    end
    
end

#######################################################################
