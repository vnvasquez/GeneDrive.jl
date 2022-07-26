using GeneDrive
using Test
using OrdinaryDiffEq
using JuMP
using Ipopt
# import LinearAlgebra

@testset "GeneDrive" begin
    include("tests_temperature_response.jl")
    include("tests_decision_model.jl")
    include("tests_dynamic_model.jl")
end
