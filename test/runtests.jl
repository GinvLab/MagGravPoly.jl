

include("mag_runtests.jl")

include("grav_runtests.jl")

include("joint_runtests.jl")


#-------------------------------------------------------#
#-------------------------TEST--------------------------#
#-------------------------------------------------------#


@testset "MagGrav2Dpoly tests" begin
    
    printstyled("\nMagnetic anomalies \n", bold=true,color=:yellow)
    
    printstyled("Forward calculation 2D \n", bold=true,color=:cyan)
    @test mag_testfwd1()
   
    printstyled("Forward calculation 2.75D \n", bold=true,color=:cyan)
    @test mag_testfwd2()

    printstyled("Mag - Gradient calculation \n", bold=true,color=:cyan)
    @test mag_testgrad()


    printstyled("\nGravity anomalies \n", bold=true,color=:yellow)

    printstyled("Forward calculation 2D \n", bold=true,color=:cyan)
    @test grav_testfwd1()
   
    printstyled("Forward calculation 2.75D \n", bold=true,color=:cyan)
    @test grav_testfwd2()

    printstyled("Gradient calculation \n", bold=true,color=:cyan)
    @test grav_testgrad()


    printstyled("\nJoint magnetic and gravity anomalies \n", bold=true,color=:yellow)

    printstyled("Forward calculation 2D \n", bold=true,color=:cyan)
    @test joint_testfwd1()
   
    printstyled("Forward calculation 2.75D \n", bold=true,color=:cyan)
    @test joint_testfwd2()

    printstyled("Gradient calculation \n", bold=true,color=:cyan)
    @test joint_testgrad()

end
