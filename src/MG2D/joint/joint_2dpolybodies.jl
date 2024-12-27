
#######################################################################

"""
$(TYPEDSIGNATURES)

 Function to return the vertical attraction & total magnetic field (2D or 2.75D) for a set of polygonal bodies defined by their corners.
Takes into account both induced and remnant magnetization.
2D formulations based on Talwani et al. (1959), Talwani & Heitzler (1964) and Blakely (1995).
2.75D formulations based on Rasmussen & Pedersen (1979) and Campbell (1983).
"""
function tjointpolybodies2D(grav_xzobs::Array{<:Real,2},mag_xzobs::Array{<:Real,2},northxax::Real,pbodies::JointPolygBodies2D)
    
    # for autodiff
    #@assert eltype(Jinds[1].mod)==eltype(bodies.geom.allvert)
    nbody = length(pbodies.geom.bo)
    tgrav = zeros(eltype(pbodies.geom.allvert),size(grav_xzobs,1))
    tmag = zeros(eltype(pbodies.geom.allvert),size(mag_xzobs,1))

    if pbodies.ylatext==nothing
        # 2D formulation
        for i=1:nbody
            tgrav .+= MagGrav2Dpoly.tgravpoly2D(grav_xzobs,pbodies.rho[i],pbodies.geom.bo[i])
            tmag .+= MagGrav2Dpoly.tmagpoly2D(mag_xzobs,pbodies.Jind.mod[i],pbodies.Jind.Ideg[i],pbodies.Jind.Ddeg[i],
                                pbodies.Jrem.mod[i],pbodies.Jrem.Ideg[i],pbodies.Jrem.Ddeg[i],
                                northxax,pbodies.geom.bo[i])
        end
        
    else
        # 2.75D formulation
        for i=1:nbody
            tgrav .+= MagGrav2Dpoly.tgravpoly2_75D(grav_xzobs,pbodies.rho[i],pbodies.geom.bo[i],pbodies.ylatext[1],pbodies.ylatext[2])
            tmag .+= MagGrav2Dpoly.tmagpoly2_75D(mag_xzobs,pbodies.Jind.mod[i],pbodies.Jind.Ideg[i],pbodies.Jind.Ddeg[i],
                                  pbodies.Jrem.mod[i],pbodies.Jrem.Ideg[i],pbodies.Jrem.Ddeg[i],
                                  northxax,pbodies.geom.bo[i],pbodies.ylatext[1],pbodies.ylatext[2])
        end
    end
    
    return tgrav,tmag
end

########################################################################################

"""
$(TYPEDSIGNATURES)

 Vertical attraction & Total magnetic field (2D or 2.75D) for a set of polygonal bodies defined by their corners.
Takes into account both induced and remnant magnetization.
Gravity calculation is based on two different algorithm formulations defined by `forwtype_grav`, passed as a string:
  - "talwani"      --> Talwani et al. (1959)
  - "wonbev"       --> Won & Bevis (1987)
Magnetic calculation instead is based on four different algorithm formulations defined by `forwtype_mag`, passed as a string:
  - "talwani"      --> Talwani & Heitzler (1964)
  - "talwani_red"  --> Talwani & Heitzler (1964) rederived from Kravchinsky et al. (2019)
  - "krav"         --> Kravchinsky et al. (2019) rectified by Ghirotto et al. (2021)
  - "wonbev"       --> Won & Bevis (1987)
"""
function tjointpolybodies2Dgen(grav_xzobs::Array{<:Real,2},mag_xzobs::Array{<:Real,2},northxax::Real,pbodies::JointPolygBodies2D,forwtype_grav::String,forwtype_mag::String)
    
    # for autodiff
    #@assert eltype(Jinds[1].mod)==eltype(bodies.geom.allvert)
    nbody = length(pbodies.geom.bo)
    tgrav = zeros(eltype(pbodies.geom.allvert),size(grav_xzobs,1))
    tmag = zeros(eltype(pbodies.geom.allvert),size(mag_xzobs,1))
    
    
    # 2D formulation based on forwtype_grav & forwtype_mag 
    for i=1:nbody
        tgrav .+= MagGrav2Dpoly.tgravpoly2Dgen(grav_xzobs,pbodies.rho[i],pbodies.geom.bo[i],forwtype_grav)
        tmag .+= MagGrav2Dpoly.tmagpoly2Dgen(mag_xzobs,pbodies.Jind.mod[i],pbodies.Jind.Ideg[i],pbodies.Jind.Ddeg[i],
                                         pbodies.Jrem.mod[i],pbodies.Jrem.Ideg[i],pbodies.Jrem.Ddeg[i],
                                         northxax,pbodies.geom.bo[i],forwtype_mag)
    end
    
    return tgrav,tmag
end

#####################################################################################################
