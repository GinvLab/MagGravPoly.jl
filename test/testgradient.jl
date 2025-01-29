
using Revise
using Mag2Dpoly
using GeoPoly
using LinearAlgebra
using ForwardDiff
using ReverseDiff
#using PyPlot

############################################################

# function comparegrads()

#     grad = computegrads()

#     figure()

#     subplot(211)
#     for k in keys(grad)
#         plot(grad[k],label=k)
#         #@show k,grad[k]
#     end
#     legend()

#     subplot(212)
#     for k in keys(grad)
#         if k!="findiff"
#             plot(grad[k]-grad["findiff"],label=string(k," - findiff"))
#         end
#     end
#     legend()

#     tight_layout()
#     savefig("comparisongradiens.pdf")
#     show()
#     return
# end

############################################################

function computegrads( )

    # create a misfit struct and a test model (derivative computed at such model)
    magmisf,refmodpar,testmodpar = defineprob()

    # gradient dictionary
    grad = Dict()

    # gradient by autodiff
    ADkindlist = ["FWDdiff","REVdiffTAPE","REVdiffTAPEcomp"]
    for ADkind in ADkindlist
        grad[ADkind] = gradautodiff(magmisf,testmodpar,ADkind)
    end
    
    # gradient by finite difference
    grad["findiff"] = gradfindiff(magmisf,testmodpar,Δh=0.001)

    return grad
end



############################################################

function gradautodiff(magmisf,vecmodpar,ADkind)

    if ADkind=="REVdiffTAPE"
        autodiffstuff = ReverseDiff.GradientTape(magmisf,vecmodpar)

    elseif ADkind=="REVdiffTAPEcomp"
        # compiled tape
        tape = ReverseDiff.GradientTape(magmisf,vecmodpar)
        autodiffstuff = ReverseDiff.compile(tape)

    elseif ADkind=="FWDdiff"
        if length(vecmodpar)<20
            chuncksize = length(vecmodpar)
        else
            chuncksize = 20 # 20
        end
        cfggrd = ForwardDiff.GradientConfig(magmisf,vecmodpar,ForwardDiff.Chunk{chuncksize}())
        autodiffstuff = cfggrd

    else
        autodiffstuff = nothing
        
    end

    gra = Mag2Dpoly.∇misf( magmisf,vecmodpar,ADkind,autodiffstuff )

    return gra
end

###############################################################

"""
 Gradient by finite differences using central differences
"""
function gradfindiff(magmisf,modpar ; Δh=0.01 )
    npar = length(modpar)
    gradFD = zeros(npar)
    
    for i=1:npar      
        #perc = round((j/grd.ny)*100,digits=2)
        #print("\r FD ",perc," % done")
        # central differences
        # 1
        curmod1 = copy(modpar)
        curmod1[i] -= Δh # -
        mv1 = magmisf(curmod1)
        # 2
        curmod2 = copy(modpar)
        curmod2[i] += Δh # +
        mv2 = magmisf(curmod2)
        # grad
        gradFD[i] = (mv2-mv1)/(2*Δh)
    end
    return gradFD
end

########################################################################

function defineprob()

    #Def. northxax, xzobs, magnetization parameters, bodyindices -> qcur, qnew
    northxax = 90.0
    N=40
    xzobs = hcat(LinRange(-100.0,300.0,N), -100.0*ones(N))

    #Def. of Topography data struct.
    topo = TopoEdges(xzobs)
    #topography[:,2] .-= 100.0 


    northxax = 90.0
    N=40
    xzobs = hcat(LinRange(-100.0,300.0,N), -100.0*ones(N))

    #Def. of Topography data struct.
    topo = TopoEdges(xzobs)
    #topography[:,2] .-= 100.0 


    ############################################################
    ##  Ref mod
    
    smallpoly = true

        
    if smallpoly

        Jindref = MagnetizVector(mod=[0.45],Ideg=[-83.0],Ddeg=[134.0])
        Jremref = MagnetizVector(mod=[3.0],Ideg=[-83.0],Ddeg=[134.0])

        verticesref=[50.0 30.0;
                     60.0 90.0;
                     80.0 90.0;
                     110.0 30.0]

        ind1 = collect(1:4)
        bodyindicesref = [ind1]

    else

        Jindref = MagnetizVector(mod=[0.01,0.45],Ideg=[-83.0,-83.0],Ddeg=[134.0,134.0])
        Jremref = MagnetizVector(mod=[5.0,3.0],Ideg=[83.0,-83.0],Ddeg=[-45.0,134.0])

        verticesref  = [35.0 50.0;
                        65.0 50.0;
                        80.0 35.0;
                        65.0 20.0;
                        35.0 20.0;
                        20.0 35.0;   #FIN QUI
                        95.0 50.0;
                        125.0 50.0;
                        140.0 35.0;
                        125.0 20.0;
                        95.0 20.0;
                        80.0 35.0]
        
        ind1 = collect(1:6)
        ind2 = collect(7:12)
        bodyindicesref = [ind1,ind2] #[ind1] 

    end

    pbodyref = MagPolygBodies2D(bodyindicesref,verticesref,Jindref,Jremref)

    verticesref,bodyindicesref = GeoPoly.checkbodyindices(verticesref,bodyindicesref)

    refmod = Mag2Dpoly.magstruct2vec(pbodyref)
    
    ############################################################
    ##  Test mod
    
    if smallpoly

        Jind = MagnetizVector(mod=[2.5],Ideg=[-63.0],Ddeg=[105.0])
        Jrem = MagnetizVector(mod=[5.0],Ideg=[-85.0],Ddeg=[110.0])  

        vertices=[40.0 40.0;
                  40.0 100.0;
                  100.0 100.0;
                  100.0 40.0]

        ind1 = collect(1:4)
        bodyindices = [ind1]

    else

        Jind = MagnetizVector(mod=[0.1,0.45],Ideg=[-83.0,-83.0],Ddeg=[134.0,134.0])
        Jrem = MagnetizVector(mod=[5.0,3.0],Ideg=[83.0,-83.0],Ddeg=[-45.0,134.0])  

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

    end    
    
    pbody = MagPolygBodies2D(bodyindices,vertices,Jind,Jrem)        

    vertices,bodyindices = checkbodyindices(vertices,bodyindices)

    GeoPoly.checkmodelizdim(topo,pbody.geom.bo,40.0)

    testmod = Mag2Dpoly.magstruct2vec(pbody)

    ############################################################
        
    tmagobs = tmagpolybodies2D(xzobs,northxax,pbodyref)

    # Calculation Data covariance matrix
    sigmad = 10.0
    Cd = Diagonal(sigmad.*ones(length(tmagobs))) 
    Cdinv = inv(Cd)

    # define mymisf structure
    magmisf = Mag2DPolyMisf(bodyindices,northxax,xzobs,tmagobs,Cdinv) 

    return magmisf,refmod,testmod
end



########################################################################

