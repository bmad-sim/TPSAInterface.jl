using TPSAInterface
using Test
using JET

@testset "TPSAInterface.jl" begin
    @testset "Code linting (JET.jl)" begin
        JET.test_package(TPSAInterface; target_defined_modules = true)
    end
    # Write your tests here.
end
