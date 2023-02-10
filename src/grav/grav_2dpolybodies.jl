
#######################################################################

"""
$(TYPEDSIGNATURES)

 Vertical attraction (2D or 2.75D) for a set of polygonal bodies defined by their corners.
2D formulation based on Talwani et al. (1959) and Blakely (1995). 
2.75D formulation based on Rasmussen & Pedersen (1979). 
"""
function tgravpolybodies2D(xzobs::Array{<:Real,2},pbodies::GravPolygBodies2D)

    nbody = length(pbodies.geom.bo)
    tgrav = zeros(eltype(pbodies.geom.allvert),size(xzobs,1))

    if pbodies.ylatext==nothing
        # 2D formulation
        for i=1:nbody
            tgrav .+= tgravpoly2D(xzobs,pbodies.rho[i],pbodies.geom.bo[i]) 
        end

    else
        # 2.75D formulation
        for i=1:nbody
            tgrav .+= tgravpoly2_75D(xzobs,pbodies.rho[i],pbodies.geom.bo[i],
                                    pbodies.ylatext[1],pbodies.ylatext[2]) 
        end

    end
    
    return tgrav
end

########################################################################################

"""
$(TYPEDSIGNATURES)

 Vertical attraction (2D) for a set of polygonal bodies defined by their corners.
Generic version containing two different algorithm formulations `forwardtype`, passed as a string:
  - "talwani"      --> Talwani et al. (1959)
  - "wonbev"       --> Won & Bevis (1987)
"""
function tgravpolybodies2Dgen(xzobs::Array{<:Real,2},pbodies::GravPolygBodies2D,
                              forwardtype::String)

    @assert body.ylatext==nothing

    # if (fieldnames(typeof(pbodies))==(:geom,:rho))||
    #     (fieldnames(typeof(pbodies))==(:geom,:Jind,:Jrem,:rho))  
    # for autodiff
    #@assert eltype(Jinds[1].mod)==eltype(bodies.geom.allvert)
    nbody = length(pbodies.geom.bo)
    tgrav = zeros(eltype(pbodies.geom.allvert),size(xzobs,1))
    for i=1:nbody
        tgrav .+= tgravpoly2Dgen(xzobs,pbodies.rho[i],pbodies.geom.bo[i],forwardtype) 
    end
    # else
    #     error("pbodies must be only of type `GravPolygBodies2D` or `JointPolygBodies2D`. Aborting")
    # end
    
    return tgrav
end

########################################################################################

"""
$(TYPEDSIGNATURES)

 Vertical attraction (2D) for a polygon defined by its corners. 
Based on Talwani et al. (1959) and Blakely (1995). 
"""
function tgravpoly2D(xzobs::Array{<:Real,2},ρ::Real,body::BodySegments2D)

    ## LOOPING on segments MUST be in ANTI-CLOCKWISE order
    aclockw = checkanticlockwiseorder(body)
    
    if !aclockw
        error("tgravpoly2D(): vertices *not* ordered anticlockwise. Aborting.")
    end
  
    ## Loop on observation points and segments
    nobs = size(xzobs,1)
    ##---------------------------------------
    # grav = Vector{eltype(body.ver1)}(undef,nobs)
    ty = [eltype(body.ver1), typeof(ρ)]
    abf = (ty.<:AbstractFloat)
    allabf = all(abf)

    if allabf
        grav = Vector{ty[1]}(undef,nobs)
    elseif !abf[1]
        grav = Vector{ty[1]}(undef,nobs)
    elseif !abf[2]
        grav = Vector{ty[2]}(undef,nobs)
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

            tsum += gravtalwani(x1,z1,x2,z2,ρ)
        end
        grav[iob] = tsum
    end
    return grav
end


########################################################################

"""
$(TYPEDSIGNATURES)

 Vertical attraction (2.5D) for a polygon defined by its corners. 
Based on Rasmussen & Pedersen (1979). 
"""
function tgravpoly2_75D(xzobs::Array{<:Real,2},ρ::Real,body::BodySegments2D,
                     y1::Real,y2::Real)

    ## LOOPING on segments MUST be in ANTI-CLOCKWISE order
    aclockw = checkanticlockwiseorder(body)
    
    if !aclockw
        error("tgravpoly2_75D(): vertices *not* ordered anticlockwise. Aborting.")
    end
  
    ## Loop on observation points and segments
    nobs = size(xzobs,1)
    ##---------------------------------------
    # grav = Vector{eltype(body.ver1)}(undef,nobs)
    ty = [eltype(body.ver1), typeof(ρ)]
    abf = (ty.<:AbstractFloat)
    allabf = all(abf)

    if allabf
        grav = Vector{ty[1]}(undef,nobs)
    elseif !abf[1]
        grav = Vector{ty[1]}(undef,nobs)
    elseif !abf[2]
        grav = Vector{ty[2]}(undef,nobs)
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

            tsum += grav2_75D(x1,y1,z1,x2,y2,z2,ρ)
        end
        grav[iob] = tsum
    end
    return grav
end


###################################################################################

"""
$(TYPEDSIGNATURES)

 Vertical attraction (2D) for a polygon defined by its corners.
Generic version containing two different algorithm formulations `forwardtype`, passed as a string:
  - "talwani"      --> Talwani et al. (1959)
  - "wonbev"       --> Won & Bevis (1987)
"""
function tgravpoly2Dgen(xzobs::Array{<:Real,2},ρ::Real,body::BodySegments2D,
                        forwardtype::String)

    ## LOOPING on segments MUST be in ANTI-CLOCKWISE order
    aclockw = checkanticlockwiseorder(body)
    
    if !aclockw
        error("tgravpoly2D(): vertices *not* ordered anticlockwise. Aborting.")
    end
  
    ## Loop on observation points and segments
    nobs = size(xzobs,1)
    ##---------------------------------------
    # grav = Vector{eltype(body.ver1)}(undef,nobs)
    ty = [eltype(body.ver1), typeof(ρ)]
    abf = (ty.<:AbstractFloat)
    allabf = all(abf)

    if allabf
        grav = Vector{ty[1]}(undef,nobs)
    elseif !abf[1]
        grav = Vector{ty[1]}(undef,nobs)
    elseif !abf[2]
        grav = Vector{ty[2]}(undef,nobs)
    end
    
    ##--------------------------------------- 
    # check right forwardtype
    if forwardtype != "talwani" && forwardtype != "wonbev"
        error("tgravpoly2Dgen(): [forwardtype] must be 'talwani' or 'wonbev'")
    end 
    
    ##-------------------------------------
    ## Loop on observation points and segments
    nobs = size(xzobs,1)
    grav = Vector{typeof(ρ)}(undef,nobs)

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
                tsum += gravtalwani(x1,z1,x2,z2,ρ)
                
            elseif forwardtype == "wonbev"
                tsum += gravwonbev(x1,z1,x2,z2,ρ)
                
            end
        end
        
        grav[iob] = tsum
    end

    return grav
end

########################################################################


"""
$(TYPEDSIGNATURES)

 Vertical attraction (2D) calculation for a polygon side defined by its corners. 
Based on Talwani et al. (1959) and Blakely (1995). 
"""
function gravtalwani(x1::Real,z1::Real,x2::Real,z2::Real,ρ::Real)

    # Quantities for errors definitions
    small = 1e4*eps(typeof(x1))
    anglelim = 0.995*π
    
    #Definition of Gravitational Constant
    γ = 6.6743e-11

    #-------------------------    
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

    denom = z2 - z1
    
    #Check on denominator ≠ 0
    if denom == 0.0
        denom = small
    end
    
    r1sq = x1^2 + z1^2
    r2sq = x2^2 + z2^2

    θdiff = atan(z2,x2) - atan(z1,x1)

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

    #--------------------------------------------
    α = (x2 - x1)/denom
    β = (x1*z2 - x2*z1)/denom
    term1 = β/(1.0 + α^2)
    term2 = 0.5*(log(r2sq) - log(r1sq))
    
    eq = term1*(term2 - α*θdiff)

    #--------------------------------------------
    #Minus sign to take into account opposite looping on polygon segments
    factor = -2.e5*ρ*γ
    
    g = factor*eq
    
    return g
end

#######################################################################################

""" 
$(TYPEDSIGNATURES)

 Vertical attraction (2D) calculation for a polygon side defined by its corners. Formulas from Won & Bevis (1987).
"""
function gravwonbev(x1::Real,z1::Real,x2::Real,z2::Real,ρ::Real)

    # Quantities for errors definitions
    small = 1e4*eps(typeof(x1))
    anglelim = 0.995*π

    #Definition of Gravitational Constant
    γ = 6.6743e-11

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
    
    R = x21^2+z21^2

    #------------------------

    r1 = x1^2+z1^2
    r2 = x2^2+z2^2

    lor21 = 0.5*log(r2) - 0.5*log(r1)
    
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
    # Calculation of g
    factor = -2.e5*ρ*γ
    
    if x21 != 0.0       
        g = factor*(x21*(x1*z2 - x2*z1)*(θdiff + 0.5*lor21*(z21/x21)))/R         
    else
        g = factor*x1*lor21      
    end
    
    return g
end
    
#######################################################################################

"""
$(TYPEDSIGNATURES)

 Vertical attraction (2.5D) calculation for a polygon side defined by its corners. 
Based on Rasmussen & Pedersen (1979). 
"""
function grav2_75D(x1::Real,y1::Real,z1::Real,x2::Real,y2::Real,z2::Real,ρ::Real)
   
    # Check on dimensions of y1 and y2
    @assert y1 < y2

    # Quantities for errors definitions
    small = 1e4*eps(typeof(x1))
    anglelim = 0.995*π
    
    #Definition of Gravitational Constant
    γ = 6.6743e-11

    #-------------------------    
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
        z21 = small
    end
    
    s = sqrt(x21^2+z21^2)

    ϕ = atan(z21,x21)
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
    lor2(y::Real) = log((u2+R2(y))/(u1+R1(y)))
    atan1(y::Real) = atan((u1*y),(w*R1(y)))
    atan2(y::Real) = atan((u2*y),(w*R2(y)))
    diffatan(y::Real) = atan2(y)-atan1(y)
    #-------------------------------

    # Calculation of terms for gravity response
    lor = lor1(y2)+lor1(-y1)
    
    diffatan2 = diffatan(y2)
    diffatan1 = diffatan(y1)

    factor = -1.e5*ρ*γ

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
    # Continuing calculation of terms for gravity response
    difatan = diffatan2 - diffatan1

    #-------------------------------------------------
    # Calculation of bodies' gravity anomaly
    factor = -1.e5*ρ*γ
    
    g = factor*(cosϕ*(y2*lor2(y2)-y1*lor2(y1))+((x1*z2-z1*x2)/s)*(cosϕ*(difatan)-sinϕ*lor))
   
    return g
end

#######################################################################################
