using GeneDrive
using Test
using OrdinaryDiffEq
# import LinearAlgebra

@testset "GeneDrive" begin
    include("tests_temperature_response.jl")
    include("tests_dynamic_model.jl")
end
