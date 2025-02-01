

########################################################################

"""
$(TYPEDSIGNATURES)

 Total magnetic field (2D or 2.75D) for a set of polygonal bodies defined by their corners. Takes into account both induced and remnant magnetization.
2D formulation based on Talwani & Heitzler (1964), the default algorithm in Mag2Dpoly package.
2.75D formulation based on Rasmussen & Pedersen (1979) and Campbell (1983).

"""
function tmagpolybodies2D(xzobs::Array{<:Real,2},northxax::Real,pbodies::MagPolygBodies2D)

    nbody = length(pbodies.geom.bo)
    tmag = zeros(eltype(pbodies.geom.allvert),size(xzobs,1))

    if pbodies.ylatext==nothing
        # 2D formulation
        for i=1:nbody
            tmag .+= tmagpoly2D(xzobs,pbodies.Jind.mod[i],pbodies.Jind.Ideg[i],pbodies.Jind.Ddeg[i],
                                pbodies.Jrem.mod[i],pbodies.Jrem.Ideg[i],pbodies.Jrem.Ddeg[i],
                                northxax,pbodies.geom.bo[i]) 
        end


    else
        # 2.75D formulation
        for i=1:nbody
            tmag .+= tmagpoly2_75D(xzobs,pbodies.Jind.mod[i],pbodies.Jind.Ideg[i],pbodies.Jind.Ddeg[i],
                                  pbodies.Jrem.mod[i],pbodies.Jrem.Ideg[i],pbodies.Jrem.Ddeg[i],
                                  northxax,pbodies.geom.bo[i],pbodies.ylatext[1],pbodies.ylatext[2]) 
        end
    end

    return tmag
end

########################################################################

"""
$(TYPEDSIGNATURES)

 Total magnetic field (2D) for a set of polygonal bodies defined by their corners. Takes into account both induced and remnant magnetization.
Generic version containing four different algorithm formulations `forwardtype`, passed as a string:
  - "talwani"      --> Talwani & Heitzler (1964)
  - "talwani_red"  --> Talwani & Heitzler (1964) rederived from Kravchinsky et al. (2019)
  - "krav"         --> Kravchinsky et al. (2019) rectified by Ghirotto et al. (2021)
  - "wonbev"       --> Won & Bevis (1987)
"""
function tmagpolybodies2Dgen(xzobs::Array{<:Real,2},northxax::Real,pbodies::MagPolygBodies2D,
                             forwardtype::String)

    @assert pbodies.ylatext==nothing

    nbody = length(pbodies.geom.bo)
    #tmag = zeros(eltype(Jinds[1].mod),size(xzobs,1))
    tmag = zeros(eltype(pbodies.geom.allvert),size(xzobs,1))
    for i=1:nbody
        tmag .+= tmagpoly2Dgen(xzobs,pbodies.Jind.mod[i],pbodies.Jind.Ideg[i],pbodies.Jind.Ddeg[i],
                               pbodies.Jrem.mod[i],pbodies.Jrem.Ideg[i],pbodies.Jrem.Ddeg[i],
                               northxax,pbodies.geom.bo[i],forwardtype)
    end
    return tmag
end


###################################################################################

"""
$(TYPEDSIGNATURES)

 Total magnetic field (2D) for a polygon defined by its corners. Takes into account both induced and remnant magnetization.
Based on Talwani & Heitzler (1964), the default algorithm in MagGrav2Dpoly package. 
"""
function tmagpoly2D(xzobs::Array{<:Real,2},Jindmod::Real,JindIdeg::Real,JindDdeg::Real,
                    Jremmod::Real,JremIdeg::Real,JremDdeg::Real,
                    northxax::Real,body::BodySegments2D)

    # for autodiff
    #@assert eltype(Jindmod)==eltype(pbodies.allvert)

    ## LOOPING on segments MUST be in ANTI-CLOCKWISE order
    aclockw = checkanticlockwiseorder(body)

    if !aclockw
        error("tmagpoly2D(): vertices *not* ordered anticlockwise. Aborting.")
    end
   
    ##---------------------
    ## Get the angles
   
    ## check modules
    @assert Jindmod >= 0.0
    @assert Jremmod >= 0.0

    ## `northxax` is the angle between geographic north and the positive x axis
    @assert 0.0 <= northxax <= 360.0
    Cnorth = deg2rad(northxax)
    
    ## check angles
    @assert -90.0 <= JindIdeg <= 90.0
    @assert -90.0 <= JremIdeg <= 90.0
    
    @assert -180.0 <= JindDdeg <= 180.0
    @assert -180.0 <= JremDdeg <= 180.0

    # deg to rad
    Iind = deg2rad(JindIdeg)
    Dind = deg2rad(JindDdeg)
    Irem = deg2rad(JremIdeg)
    Drem = deg2rad(JremDdeg)

    # Calculation of Jx and Jz only for case != wonbev
    Jtotx,_,Jtotz = magcomp(Jindmod,Iind,Dind,Jremmod,Irem,Drem,Cnorth)

    ## Loop on observation points and segments
    nobs = size(xzobs,1)
    ##---------------------------------------
    # totfield = Vector{eltype(body.ver1)}(undef,nobs)
    ty = [eltype(body.ver1), typeof(Jindmod)]
    abf = (ty.<:AbstractFloat)
    allabf = all(abf)

    if allabf
        totfield = Vector{ty[1]}(undef,nobs)
    elseif !abf[1]
        totfield = Vector{ty[1]}(undef,nobs)
    elseif !abf[2]
        totfield = Vector{ty[2]}(undef,nobs)
    end
    ##---------------------------------------

    ## loop on observation points
    for iob=1:nobs      
        
        xo = xzobs[iob,1]
        zo = xzobs[iob,2]

        ## loop on segments
        tsum = 0.0
        for ise=1:body.nsegm

            x1 = body.ver1[ise,1]-xo
            z1 = body.ver1[ise,2]-zo
            x2 = body.ver2[ise,1]-xo
            z2 = body.ver2[ise,2]-zo

            tsum += tmagtalwani(x1,z1,x2,z2,Jtotx,Jtotz,Iind,Dind,Cnorth)
        end

        # #---------------------------------------
        # # !!TESTING!!
        # x1 = body.ver1[:,1] .- xo
        # z1 = body.ver1[:,2] .- zo
        # x2 = body.ver2[:,1] .- xo
        # z2 = body.ver2[:,2] .- zo
        # tall = tmagtalwani.(x1,z1,x2,z2,Jtotx,Jtotz,Iind,Dind,Cnorth)
        # tsum2 = sum(tall)
        # #@assert tsum2≈tsum
        # #---------------------------------------
        
        # convert field from A/m to nT
        totfield[iob] = convert_H_to_B_nT( tsum )
    end

    return totfield
end

#####################################################################################

"""
$(TYPEDSIGNATURES)

 Total magnetic field (2.75D) for a polygon defined by its corners. Takes into account both induced and remnant magnetizations.
Based on Rasmussen & Pedersen (1979) and Campbell (1983). 
"""
function tmagpoly2_75D(xzobs::Array{<:Real,2},Jindmod::Real,JindIdeg::Real,JindDdeg::Real,
                      Jremmod::Real,JremIdeg::Real,JremDdeg::Real,
                      northxax::Real,body::BodySegments2D,y1::Real,y2::Real)

    # for autodiff
    #@assert eltype(Jindmod)==eltype(pbodies.allvert)

    ## LOOPING on segments MUST be in ANTI-CLOCKWISE order
    aclockw = checkanticlockwiseorder(body)

    if !aclockw
        error("tmagpoly2_75D(): vertices *not* ordered anticlockwise. Aborting.")
    end
   
    ##---------------------
    ## Get the angles
   
    ## check modules
    @assert Jindmod >= 0.0
    @assert Jremmod >= 0.0

    ## `northxax` is the angle between geographic north and the positive x axis
    @assert 0.0 <= northxax <= 360.0
    Cnorth = deg2rad(northxax)
    
    ## check angles
    @assert -90.0 <= JindIdeg <= 90.0
    @assert -90.0 <= JremIdeg <= 90.0
    
    @assert -180.0 <= JindDdeg <= 180.0
    @assert -180.0 <= JremDdeg <= 180.0

    # deg to rad
    Iind = deg2rad(JindIdeg)
    Dind = deg2rad(JindDdeg)
    Irem = deg2rad(JremIdeg)
    Drem = deg2rad(JremDdeg)

    # Calculation of Jx and Jz only for case != wonbev
    Jtotx,Jtoty,Jtotz = magcomp(Jindmod,Iind,Dind,Jremmod,Irem,Drem,Cnorth)

    ## Loop on observation points and segments
    nobs = size(xzobs,1)
    ##---------------------------------------
    # totfield = Vector{eltype(body.ver1)}(undef,nobs)
    ty = [eltype(body.ver1), typeof(Jindmod)]
    abf = (ty.<:AbstractFloat)
    allabf = all(abf)

    if allabf
        totfield = Vector{ty[1]}(undef,nobs)
    elseif !abf[1]
        totfield = Vector{ty[1]}(undef,nobs)
    elseif !abf[2]
        totfield = Vector{ty[2]}(undef,nobs)
    end
    ##---------------------------------------

    ## loop on observation points
    for iob=1:nobs      
        
        xo = xzobs[iob,1]
        zo = xzobs[iob,2]

        ## loop on segments
        tsum = 0.0
        for ise=1:body.nsegm

            x1 = body.ver1[ise,1]-xo
            z1 = body.ver1[ise,2]-zo
            x2 = body.ver2[ise,1]-xo
            z2 = body.ver2[ise,2]-zo

            tsum += tmag2_75D(x1,y1,z1,x2,y2,z2,Jtotx,Jtoty,Jtotz,Iind,Dind,Cnorth)
        end

        # #---------------------------------------
        # # !!TESTING!!
        # x1 = body.ver1[:,1] .- xo
        # z1 = body.ver1[:,2] .- zo
        # x2 = body.ver2[:,1] .- xo
        # z2 = body.ver2[:,2] .- zo
        # tall = tmagtalwani.(x1,z1,x2,z2,Jtotx,Jtotz,Iind,Dind,Cnorth)
        # tsum2 = sum(tall)
        # #@assert tsum2≈tsum
        # #---------------------------------------
        
        # convert field from A/m to nT
        totfield[iob] = convert_H_to_B_nT( tsum )
    end

    return totfield
end

#####################################################################################

"""
$(TYPEDSIGNATURES)

 Total magnetic field (2D) calculation for a polygon defined by its corners. Takes into account both induced and remnant magnetizations.
Generic version containing four different algorithm formulations `forwardtype`, passed as a string:
  - "talwani"      --> Talwani & Heitzler (1964)
  - "talwani_red"  --> Talwani & Heitzler (1964) rederived from Kravchinsky et al. (2019)
  - "krav"         --> Kravchinsky et al. (2019) rectified by Ghirotto et al. (2021)
  - "wonbev"       --> Won & Bevis (1987)
"""
function tmagpoly2Dgen(xzobs::Array{<:Real,2},Jindmod::Real,JindIdeg::Real,JindDdeg::Real,Jremmod::Real,JremIdeg::Real,JremDdeg::Real,
                       northxax::Real,body::BodySegments2D,forwardtype::String)

    # for autodiff
    #@assert eltype(Jind.mod)==eltype(pbodies.allvert)

    ## LOOPING on segments MUST be in ANTI-CLOCKWISE order
    aclockw = checkanticlockwiseorder(body)

    if !aclockw
        error("tmagpoly2Dgen(): vertices *not* ordered anticlockwise. Aborting.")
    end
   
    ##---------------------
    ## Get the angles
   
    ## check modules
    @assert Jindmod >= 0.0
    @assert Jremmod >= 0.0

    ## `northxax` is the angle between geographic north and the positive x axis
    @assert 0.0 <= northxax <= 360.0
    Cnorth = deg2rad(northxax)
    
    ## check angles
    @assert -90.0 <= JindIdeg <= 90.0
    @assert -90.0 <= JremIdeg <= 90.0
    
    @assert -180.0 <= JindDdeg <= 180.0
    @assert -180.0 <= JremDdeg <= 180.0

    # check right forwardtype
    if forwardtype != "talwani" && forwardtype != "talwani_red" && forwardtype != "krav" && forwardtype != "wonbev"
        error("tmagpoly2Dgen(): [forwardtype] must be 'talwani' or 'talwani_red' or 'krav' or 'wonbev'")
    end 
    
    # deg to rad
    Iind = deg2rad(JindIdeg)
    Dind = deg2rad(JindDdeg)
    Irem = deg2rad(JremIdeg)
    Drem = deg2rad(JremDdeg)

    # Calculation of Jx and Jz only for case != wonbev
    if forwardtype != "wonbev"
        Jtotx,_,Jtotz = magcomp(Jindmod,Iind,Dind,Jremmod,Irem,Drem,Cnorth)
    end
    
    ##-------------------------------------
    ## Loop on observation points and segments
    nobs = size(xzobs,1)
    totfield = Vector{typeof(Jindmod)}(undef,nobs)

    ## loop on observation points
    for iob=1:nobs      
        
        xo = xzobs[iob,1]
        zo = xzobs[iob,2]

        ## loop on segments
        tsum = 0.0
        for ise=1:body.nsegm

            x1 = body.ver1[ise,1]-xo
            z1 = body.ver1[ise,2]-zo
            x2 = body.ver2[ise,1]-xo
            z2 = body.ver2[ise,2]-zo

            if forwardtype == "talwani"
                tsum += tmagtalwani(x1,z1,x2,z2,Jtotx,Jtotz,Iind,Dind,Cnorth)

            elseif forwardtype == "talwani_red"
                tsum += tmagtalwanired(x1,z1,x2,z2,Jtotx,Jtotz,Iind,Dind,Cnorth)

            elseif forwardtype == "krav"
                tsum += tmagkrav(x1,z1,x2,z2,Jtotx,Jtotz,Iind,Dind,Cnorth)

            elseif forwardtype == "wonbev"
                tsum += tmagwonbev(x1,z1,x2,z2,Jindmod,Jremmod,Iind,Dind,Irem,Drem,Cnorth)
                
            end            
         end                   
                          
        # convert field from A/m to nT
        totfield[iob] = convert_H_to_B_nT( tsum )
    end

    return totfield
end

########################################################################

"""
$(TYPEDSIGNATURES)
 
 Total magnetic field (2D) calculation for a polygon side defined by its corners. Takes into account both induced and remnant magnetizations. Formulas from Talwani & Heitzler (1964).
"""
function tmagtalwani(x1::Real,z1::Real,x2::Real,z2::Real,
                     Jx::Real,Jz::Real,Iind::Real,Dind::Real,C::Real)
    
    # Quantities for error definitions
    small = 1e4*eps(typeof(x1))
    anglelim = 0.995*π

    #----------------------------------
    # Error if a corner is too close to the observation point (calculation continues)
    # and the corner are slightly moved away        
    if abs(x1) < small && abs(z1) < small
        x1 = flipsign(small,x1)
        z1 = flipsign(small,z1)
        @warn "A corner is too close to an observation point (calculation continues)"
    end
    
    if abs(x2) < small && abs(z2) < small
        x2 = flipsign(small,x2)
        z2 = flipsign(small,z2)
        @warn "A corner is too close to an observation point (calculation continues)"
    end
    
    #--------------
    x21 = x2-x1
    z21 = z2-z1
    
    # Check on z21 ≠ 0
    if z21 == 0.0
        ϕ = -acot(x21/small)
    else
        ϕ = -acot(x21/z21)
    end
       
    #--------------------------------------    
    r1 = x1^2 + z1^2
    r2 = x2^2 + z2^2
    
    flog = 0.5*(log(r2)-log(r1))
        
    #-----------------------
    # Get the angles
    θ1 = atan(z1,x1)
    θ2 = atan(z2,x2)
    θdiff = θ2-θ1

    
    # In the case polygon sides cross the x axis
    if θdiff < -π
        θdiff = θdiff + 2.0*π
    elseif θdiff > π
        θdiff = θdiff - 2.0*π
    end   
    
    
    # Error if the side is too close to the observation point (calculation continues)
    if abs(θdiff) > anglelim
        @warn "A polygon side is too close to an observation point (calculation continues)"
    end
    
    #----------------------------------------------------- 
    
    # vertical component
    V = 2.0*sin(ϕ) * (Jx*( (θdiff)*cos(ϕ) + sin(ϕ)*flog) -
                      Jz*( (θdiff)*sin(ϕ) - cos(ϕ)*flog) )

    # horizonatal component
    H = 2.0*sin(ϕ) * (Jx*( (θdiff)*sin(ϕ) - cos(ϕ)*flog) +
                      Jz*( (θdiff)*cos(ϕ) + sin(ϕ)*flog) )

    #-------------------------------------------------
    # Divided by 4π to obtain a magnetic anomaly value in nT (magnetization given in A/m) 
    totfield = (1.0/(4.0*π)) * (H*cos(Iind)*cos(C-Dind) + V*sin(Iind))

    return totfield 
end

#####################################################################################

"""
$(TYPEDSIGNATURES)
 
 Total magnetic field (2.75D) calculation for a polygon side defined by its corners. Takes into account both induced and remnant magnetizations. Formulas derived from Rasmussen & Pedersen (1979) and Campbell (1983).
"""
function tmag2_75D(x1::Real,y1::Real,z1::Real,x2::Real,y2::Real,z2::Real,
                  Jx::Real,Jy::Real,Jz::Real,Iind::Real,Dind::Real,C::Real)

    # Check on dimensions of y1 and y2
    @assert y1 < y2
    
    # Quantities for error definitions
    small = 1e4*eps(typeof(x1))
    anglelim = 0.995*π
    
    #--------------------------------------
    # Error if a corner is too close to the observation point (calculation continues)
    # and the corner are slightly moved away
    if abs(x1) < small && abs(z1) < small
        x1 = flipsign(small,x1)
        z1 = flipsign(small,z1)
        @warn "A corner is too close to an observation point (calculation continues)"
    end
    
    if abs(x2) < small && abs(z2) < small
        x2 = flipsign(small,x2)
        z2 = flipsign(small,z2)
        @warn "A corner is too close to an observation point (calculation continues)"
    end
    
    #-------------------------------------   
    ## Definition of some quantities    
    x21 = x2-x1
    z21 = z2-z1
    
    # Check on z21 ≠ 0
    if z21 == 0.0
        ϕ = atan(small,x21)
    else
        ϕ = atan(z21,x21)
    end
    
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ) 
    
    u1 = cosϕ*x1 + sinϕ*z1 
    u2 = cosϕ*x2 + sinϕ*z2
    w = -sinϕ*x1 + cosϕ*z1
    
    r1 = sqrt(u1^2+w^2)
    r2 = sqrt(u2^2+w^2)

    #-------------------------------
    # Definition of some functions
    R1(y::Real) = sqrt(u1^2+w^2+y^2)
    R2(y::Real) = sqrt(u2^2+w^2+y^2)
    lor1(y::Real) = log((r1*(R2(y)+y))/(r2*(R1(y)+y)))
    lor2(y::Real) = log(((R2(y)-u2)*(R1(y)+u1))/((R2(y)+u2)*(R1(y)-u1)))
    atan1(y::Real) = atan((u1*y),(w*R1(y)))
    atan2(y::Real) = atan((u2*y),(w*R2(y)))
    diffatan(y::Real) = atan2(y)-atan1(y)
    #-------------------------------
    
    # Calculation of terms for magnetic response
    lor = lor1(y2) + lor1(-y1)

    diffatan2 = diffatan(y2)
    diffatan1 = diffatan(y1)

    #-------------------------------------------------
    # In the case polygon sides cross the x axis    
    if diffatan2 < -π
        diffatan2 = diffatan2 + 2.0*π
    elseif diffatan2 > π
        diffatan2 = diffatan2 - 2.0*π
    end    

    if diffatan1 < -π
        diffatan1 = diffatan1 + 2.0*π
    elseif diffatan1 > π
        diffatan1 = diffatan1 - 2.0*π
    end
    
    # Error if the side is too close to the observation point (calculation continues)
    if abs(diffatan2) > anglelim && abs(diffatan1) > anglelim 
        @warn "A polygon side is too close to an observation point (calculation continues)"
    end
    
    #-------------------------------------------------
    # Continuing calculation of terms for magnetic response
    difatan = diffatan2 - diffatan1
    
    term1 = cosϕ*Jx + sinϕ*Jz
    term2 = cosϕ*Jz - sinϕ*Jx
    term3 = 0.5*lor2(y1) - 0.5*lor2(y2)
    term4 = term1*lor - term2*difatan - Jy*term3
      
    #-------------------------------------------------
    # Calculation of bodies' magnetic anomaly components
    # x component
    Bx = -sinϕ*term4  
    # z component
    Bz = cosϕ*term4
    # y component
    By = (sinϕ*Jx - cosϕ*Jz)*term3 + Jy*difatan
    
    #-------------------------------------------------
    # Divided by 4π to obtain a magnetic anomaly value in nT (magnetization given in A/m)
    totfield = (1.0/(4.0*π)) * (Bx*cos(Iind)*cos(C-Dind) + Bz*sin(Iind) + By*cos(Iind)*sin(C-Dind))
    
    # Minus sign to take into account opposite looping on polygon segments
    return totfield 
end

###########################################################################

"""
$(TYPEDSIGNATURES)

 Total magnetic field (2D) calculation for a polygon side defined by its corners. Takes into account both induced and remnant magnetizations. Formulas from Kravchinsky et al (2019) rectified by Ghirotto et al. (2021). 
"""
function tmagkrav(x1::Real,z1::Real,x2::Real,z2::Real,
                  Jtotx::Real,Jtotz::Real,Iind::Real,Dind::Real,Cnorth::Real)


    # Quantities for errors definitions
    small = 1e4*eps(typeof(x1))
    anglelim = 0.995*π

    #----------------------------------
    # Error if a corner is too close to the observation point (calculation continues)
    # and the corner are slightly moved away
    if abs(x1) < small && abs(z1) < small
        x1 = flipsign(small,x1)
        z1 = flipsign(small,z1)
        @warn "A corner is too close to an observation point (calculation continues)"
    end
    
    if abs(x2) < small && abs(z2) < small
        x2 = flipsign(small,x2)
        z2 = flipsign(small,z2)
        @warn "A corner is too close to an observation point (calculation continues)"
    end
    
    #--------------
    x21 = x2-x1
    z21 = z2-z1

    #Check on z21 ≠ 0
    if z21 == 0.0
        z21 = small
    end
    
    #-------------
    tmpγ = sqrt(x21^2+z21^2)
    γx = x21 / tmpγ
    γz = z21 / tmpγ

    g = x21/z21

    if x1 >= g*z1  
        δ = 1.0
    elseif x1 < g*z1
        δ = -1.0
    end

    #--------------------
    r1 = x1^2 + z1^2
    r2 = x2^2 + z2^2
   
    lor21 = 0.5*(log(r2)-log(r1))
    
    #--------------------
    # Get the angles
    α1 = atan(δ*(z1+g*x1),(x1-g*z1))
    α2 = atan(δ*(z2+g*x2),(x2-g*z2))
    
    
    #In the case polygon sides cross the x axis
    αdiff = α2 - α1
    
    if αdiff < -π
        αdiff = αdiff + 2.0*π
    elseif αdiff > π
        αdiff = αdiff - 2.0*π
    end
    
    #----------------------------------
    
    # Error if the side is too close to the observation point (calculation continues)
    if abs(αdiff) > anglelim
        @warn "A polygon side is too close to an observation point (calculation continues)"
    end
           
    #------------------------------------
    
    P = γz*γx*lor21 + δ*(γz^2)*(αdiff)
    Q = (γz^2)*lor21 - δ*γx*γz*(αdiff)
    
    ## horizonatl and vertical field components
    H = 1.0/(2.0*pi) * (Jtotz*Q + Jtotx*P)
    V = 1.0/(2.0*pi) * (Jtotx*Q - Jtotz*P)

    #-------------------------------------------------
    ## total field anomaly 
    totfield = V*sin(Iind)+H*cos(Iind)*cos(Cnorth-Dind)
    
    return totfield
end

########################################################################################

""" 
$(TYPEDSIGNATURES)

 Total magnetic field (2D) calculation for a polygon side defined by its corners. Takes into account both induced and remnant magnetizations. Formulas from Talwani & Heitzler (1964) modified by Kravchinsky et al. (2019).
"""
function tmagtalwanired(x1::Real,z1::Real,x2::Real,z2::Real,
                        Jx::Real,Jz::Real,Iind::Real,Dind::Real,C::Real)


    # Quantities for errors definitions
    small = 1e4*eps(typeof(x1))
    anglelim = 0.995*π

    #----------------------------------
    # Error if a corner is too close to the observation point (calculation continues)
    # and the corner are slightly moved away
    if abs(x1) < small && abs(z1) < small
        x1 = flipsign(small,x1)
        z1 = flipsign(small,z1)
        @warn "A corner is too close to an observation point (calculation continues)"
    end
    
    if abs(x2) < small && abs(z2) < small
        x2 = flipsign(small,x2)
        z2 = flipsign(small,z2)
        @warn "A corner is too close to an observation point (calculation continues)"
    end
    
    #--------------
    x21 = x2-x1
    z21 = z2-z1

    # Check on z21 ≠ 0
    if z21 == 0.0    
        ϕ = -acot(x21/small)
    else
        ϕ = -acot(x21/z21)
    end
    
    #-----------------------   
    
    r1 = x1^2 + z1^2
    r2 = x2^2 + z2^2

    flog = 0.5*(log(r2)-log(r1))

    #---------------------
    
    den1 = x1+z1*cot(ϕ)
    den2 = x2+z2*cot(ϕ)
    num1 = z1-x1*cot(ϕ)
    num2 = z2-x2*cot(ϕ)

    # Controls on signs of atan argument (abs in den1 and den2)
    #-----------------------
    if den1 < 0.0
        den1 = -den1
        δ = -1.0
        θ1 = atan(num1,den1)
    else
        δ = 1.0
        θ1 = atan(num1,den1)
    end 

    if den2 < 0.0
        den2 = -den2
        θ2 = atan(num2,den2)
    else
        θ2 = atan(num2,den2)        
    end
    #-----------------------

    # In the case polygon sides cross the x axis
    θdiff = θ2-θ1
    
    if θdiff < -π
        θdiff = θdiff + 2.0*π
    elseif θdiff > π
        θdiff = θdiff - 2.0*π
    end    
    
    #-------------------------------------------------
    
    # Error if the side is too close to the observation point (calculation continues)
    if abs(θdiff) > anglelim
        @warn "A polygon side is too close to an observation point (calculation continues)"
    end
    
    #----------------------------------------------------

    # vertical component    
    V = 2.0*sin(ϕ) * (Jx * (δ*(θdiff)*cos(ϕ) + sin(ϕ)*flog)-
                      Jz * (δ*(θdiff)*sin(ϕ) - cos(ϕ)*flog) )
 
    # horizontal component
    H = 2.0*sin(ϕ) * (Jx * (δ*(θdiff)*sin(ϕ) - cos(ϕ)*flog)+
                      Jz * (δ*(θdiff)*cos(ϕ) + sin(ϕ)*flog) )

    #-------------------------------------------------
    # Divided by 4π to obtain a magnetic anomaly value in nT (magnetization given in A/m)
    totfield = (1.0/(4.0*π)) * (H*cos(Iind)*cos(C-Dind) + V*sin(Iind))
    
    return totfield  
end


#######################################################################################################################

""" 
$(TYPEDSIGNATURES)

 Total magnetic field (2D) calculation for a polygon side defined by its corners. Takes into account both induced and remnant magnetizations. Formulas from Won & Bevis (1987).
"""
function tmagwonbev(x1::Real,z1::Real,x2::Real,z2::Real,
                    modJind::Real,modJrem::Real,Iind::Real,Dind::Real,
                    Irem::Real,Drem::Real,C::Real)

    # β is angle among North and profle direction
    βi = Dind - C + π/2
    βr = Drem - C + π/2
    
    #-------------------

    # Quantities for errors definitions
    small = 1e4*eps(typeof(x1))
    anglelim = 0.995*π
    
    #----------------------------------
    # Error if a corner is too close to the observation point (calculation continues)
    # and the corner are slightly moved away
    if abs(x1) < small && abs(z1) < small
        x1 = flipsign(small,x1)
        z1 = flipsign(small,z1)
        @warn "A corner is too close to an observation point (calculation continues)"
    end
    
    if abs(x2) < small && abs(z2) < small
        x2 = flipsign(small,x2)
        z2 = flipsign(small,z2)
        @warn "A corner is too close to an observation point (calculation continues)"
    end

    #----------------------
    x21 = x2-x1
    z21 = z2-z1
        
    # Check on z21 ≠ 0
    if z21 == 0.0
        z21 = small
    end
    
    R = x21^2 + z21^2

    #------------------------

    r1 = x1^2 + z1^2
    r2 = x2^2 + z2^2

    lor21 = 0.5*(log(r2)-log(r1))

    #------------------------
    θ1 = atan(z1,x1) 
    θ2 = atan(z2,x2)
    
    # In the case polygon sides cross the x axis
    if sign(z1) != sign(z2)
        test = x1*z2 - x2*z1
        if test > 0.0
            if z1 >= 0.0 
                θ2 = θ2 + 2.0*π
            end
        elseif test < 0.0
            if z2 >= 0.0
                θ1 = θ1 + 2.0*π
            end
        end
    end   
    
    
    # Error if the side is too close to the observation point (calculation continues)
    θdiff = θ1-θ2
    if abs(θdiff) > anglelim
        @warn "A polygon side is too close to an observation point (calculation continues)"
    end
    
    #-----------------------
    
    P = (1/R)*(x1*z2 - x2*z1)*(((x1*x21 - z1*z21)/r1)-
                                 ((x2*x21 - z2*z21)/r2))

    Q = (1/R)*(x1*z2 - x2*z1)*(((x1*z21 + z1*x21)/r1)-
                                 ((x2*z21 + z2*x21)/r2))
                                      
    if x21 != 0.0
        
        g = z21/x21
        derZz = ((x21^2)/R)*(θdiff + g*lor21) - P
        derZx = -((x21*z21)/R)*(θdiff + g*lor21) + Q
        derXz = -((x21^2)/R)*(g*θdiff - lor21) + Q
        derXx = ((x21*z21)/R)*(g*θdiff - lor21) + P
    
    else

        derZz = -P
        derZx = -((z21^2)/R)*lor21 + Q
        derXz = Q
        derXx = ((z21^2)/R)*θdiff + P
        
    end

    # Magnetic strenght components due to induced magnetization
    ΔHzind = -2.0*modJind*(sin(Iind)*derZz + sin(βi)*cos(Iind)*derZx) 
    ΔHxind = -2.0*modJind*(sin(Iind)*derXz + sin(βi)*cos(Iind)*derXx) 

    # Magnetic strenght components due to remnant magnetization
    ΔHzrem = -2.0*modJrem*(sin(Irem)*derZz + sin(βr)*cos(Irem)*derZx) 
    ΔHxrem = -2.0*modJrem*(sin(Irem)*derXz + sin(βr)*cos(Irem)*derXx) 

    ΔHztot = ΔHzind + ΔHzrem
    ΔHxtot = ΔHxind + ΔHxrem

    #-------------------------------------------------
    # Divided by 4π to obtain a magnetic anomaly value in nT (magnetization given in A/m)
    ΔHtot = (1.0/(4.0*π))*(ΔHztot*sin(Iind) + ΔHxtot*sin(βi)*cos(Iind))

    return ΔHtot
end

###################################################################################
