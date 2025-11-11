
#--------------------

using Revise
using MagGravPoly
using LinearAlgebra

#--------------------

function computegradsdcshift( )

    # create a misfit struct and a test model (derivative computed at such model)
    magmisf,testmodpar = defineprob()

    # gradient dictionary
    grad = Dict()

    #=
    # gradient by autodiff
    ADkindlist = ["FWDdiff","REVdiffTAPE","REVdiffTAPEcomp"]
    for ADkind in ADkindlist
        grad[ADkind] = gradautodiff(magmisf,testmodpar,ADkind)
    end
    =#

    #-- ∇(dcshift) --# 
    # gradient by finite difference
    grad["∇_dcshiftDF"] = MG2D.magFDdcshift(magmisf,testmodpar,Δh=0.001)

    # gradient by analitic formulae
    grad["∇_dcshiftAN"] = MG2D.magANdcshift(magmisf,testmodpar,Δh=0.001)

    #-- ∇(misfit) --#
    autodiffstuff = MG2D.precalcADstuffmag(magmisf,"FWDdiff",testmodpar)
    grad["∇_misfit"] = MG2D.calc∇misfmag(magmisf,testmodpar,"FWDdiff",autodiffstuff) 
    
    return grad
end



############################################################


function defineprob()

    #Def. northxax, xzobs, magnetization parameters, bodyindices -> qcur, qnew
    northxax = 90.0
    N=40
    xzobs = hcat(LinRange(-100.0,300.0,N), -100.0*ones(N))

    #Def. of Topography data struct.
    topo = MagGravPoly.GeoPoly.TopoEdges(xzobs)
    #topography[:,2] .-= 100.0 

    
    ############################################################
    ##  Test mod

    Jind = MG2D.MagnetizVector(mod=[0.1,0.45],Ideg=[-83.0,-83.0],Ddeg=[134.0,134.0])
    Jrem = MG2D.MagnetizVector(mod=[5.0,3.0],Ideg=[83.0,-83.0],Ddeg=[-45.0,134.0])  
    
    vertices  = [55.0 50.0;
                 85.0 50.0;
                 100.0 35.0;
                 85.0 20.0;
                 55.0 20.0;
                 40.0 35.0;  
                 85.0 50.0; 
                 145.0 50.0;
                 160.0 35.0;
                 145.0 20.0;
                 115.0 20.0;
                 100.0 35.0]

    ind1 = collect(1:6)
    ind2 = collect(7:12)
    bodyindices = [ind1,ind2] #[ind1]
    
      
    bodyindices,vertices = MagGravPoly.GeoPoly.checkbodyindices(bodyindices,vertices)
    pbody = MG2D.MagPolygBodies2D(bodyindices,vertices,Jind,Jrem,ylatext=nothing)        
    
    MagGravPoly.GeoPoly.checkmodelizdim(topo,pbody.geom.bo,40.0)

    testmod = MG2D.magstruct2vec(:vertices,pbody)
    
    ############################################################
        
    tmagobs = MG2D.tmagpolybodies2D(xzobs,northxax,pbody).+(10.0.*ones(N)) #randn(N)

    # Calculation Data covariance matrix
    sigmad = 10.0
    Cd = Diagonal(sigmad.*ones(length(tmagobs))) 
    Cdinv = inv(Cd)

    # define mymisf structure
    magmisf = MG2D.Mag2DPolyMisf(bodyindices,northxax,xzobs,tmagobs,Cdinv,:vertices,Jind=Jind,Jrem=Jrem,typemisf=:dcshift) 

    return magmisf,testmod
end



########################################################################

