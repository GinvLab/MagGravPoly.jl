

include("poly_runtests.jl")


#-------------------------------------------------------#
#-------------------------TEST--------------------------#
#-------------------------------------------------------#


@testset "GeoPolygons tests" begin
    
    printstyled("\nTesting polygon assembly \n", bold=true,color=:yellow)
    
    printstyled("Single polygon \n", bold=true,color=:cyan)
    @test poly_assembly1()

    printstyled("Two polygons \n", bold=true,color=:cyan)
    @test poly_assembly2()

end
