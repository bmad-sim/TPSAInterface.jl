using TPSAInterface
using Test
using JET

@testset "TPSAInterface.jl" begin
    @testset "Code linting (JET.jl)" begin
        JET.test_package(TPSAInterface; targetdefined_modules = true)
    end
    # Write your tests here.
end
