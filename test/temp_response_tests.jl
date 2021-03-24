@testset "Temp Response Functions" begin
    test_example_reponse = ExampleResponse(10.0, 1.0)
    response = get_temperature_response(27.0, test_example_reponse)
    @test isapprox(response, 999.0; atol = 1e-3)
end