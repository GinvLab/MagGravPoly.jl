#####################################################################################################
"""
$(TYPEDSIGNATURES)

Function to check if there is intersection beetween two segments, defined by their vertices.
The output returns the (x,z) coordinates of the point in the case of intersection.
"""
function Inter2Segm(P1x::Float64,P1z::Float64,P2x::Float64,P2z::Float64,
                    P3x::Float64,P3z::Float64,P4x::Float64,P4z::Float64)


    ##Metodo di risoluzione mediante Cramer

    #Calcolo Delta0
    Delta0 = (P4z - P3z) * (P2x - P1x) - (P4x - P3x) * (P2z - P1z)
    #Se delta0 e' nullo allora trattasi di rette parallele
    if Delta0 == 0.0
        return nothing
        ## Vertici coincidenti
    elseif (P1x == P3x && P1z == P3z) ||
        (P1x == P4x && P1z == P4z) ||
        (P2x == P3x && P2z == P3z) ||
        (P2x == P4x && P2z == P4z)     
        return nothing
    else
        Delta1 = (P4x - P3x) * (P1z - P3z) - (P4z - P3z) * (P1x - P3x)
        Delta2 = (P2x - P1x) * (P1z - P3z) - (P2z - P1z) * (P1x - P3x)
        ka = Delta1 / Delta0
        kb = Delta2 / Delta0
        
        if (ka >= 0 && ka <= 1) && (kb >= 0 && kb <= 1)
            P0x = round((P1x + ka * (P2x - P1x)),digits=3)
            P0z = round((P1z + ka * (P2z - P1z)),digits=3)
            
            ## Estremità segmento cade all'interno dell'altro segmento
            if isapprox([P0x,P0z],[P1x,P1z],atol=1e-3) ||
                isapprox([P0x,P0z],[P2x,P2z],atol=1e-3) ||
                isapprox([P0x,P0z],[P3x,P3z],atol=1e-3) ||
                isapprox([P0x,P0z],[P4x,P4z],atol=1e-3)
            return nothing
            end            
        else
            return nothing
        end
    end
    P0 = hcat(P0x,P0z)
    return P0
           
end #function

###############################################################

"""
$(TYPEDSIGNATURES)

Function to check if a point is internal to a polygon.
"""
function isInternal(bo::BodySegments2D,px::Float64,pz::Float64)
    # ray-casting algorithm based on
    # https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html
    
    inside = false
    
    #Check if the point belongs to the body vertices
    for l=1:length(bo.ver1[:,1])
        #if (([px,pz].== bo.ver1[l,:]) == ones(2)) == true
        if all([px,pz].== bo.ver1[l,:])
            return inside
        end
    end
    
    j=length(bo.ver1[:,1])
    
    for i=1:length(bo.ver1[:,1])    

        xi = bo.ver1[i,1]
        zi = bo.ver1[i,2]
        xj = bo.ver1[j,1]
        zj = bo.ver1[j,2]
        
        intersect = ((zi>pz)!=(zj>pz))&&(px<(xj-xi)*(pz-zi)/(zj-zi)+xi)
        if intersect
            inside = !inside
        end
        j=i
    end 
    return inside
end

###############################################################

