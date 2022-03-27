@testset "Temp Response Functions" begin

    # NoTemp
    test_notemp = NoResponse()
    response = get_temperature_response(0.0, test_notemp)
    @test response == 0.0

    # AedesAegypti: Rossi et al 2014
    test_eggmortalityrossi = EggMortalityRossi(0.0731, 0.0595)
    response = get_temperature_response(27.0, test_eggmortalityrossi)
    @test isapprox(response, 0.364427; atol=1e-3)

    test_eggdurationrossi =
        EggDurationRossi(0.00764, 273.0, 40.55, 13094.10, 92.501, 28169.2)
    response = get_temperature_response(27.0, test_eggdurationrossi)
    @test isapprox(response, 0.670839; atol=1e-3)

    #= #TODO: remove Pelotti implementation once decided
    test_eggdurationrossi = EggDurationRossi(6.9, 4.0, 20.0, 4.1)
    response = get_temperature_response(27.0, test_eggdurationrossi)
    @test isapprox(response, 0.1496296; atol = 1e-3)
    =#

    test_larvamortalityrossi = LarvaMortalityRossi(0.0143, 0.00189)
    response = get_temperature_response(27.0, test_larvamortalityrossi)
    @test isapprox(response, 0.056716; atol=1e-3)

    test_larvadurationrossi = LarvaDurationRossi(0.00219, 273.0, 25.21, 7514.34)
    response = get_temperature_response(27.0, test_larvadurationrossi)
    @test isapprox(response, 0.772694; atol=1e-3)

    test_pupamortalityrossi = PupaMortalityRossi(0.0143, 0.00189)
    response = get_temperature_response(27.0, test_pupamortalityrossi)
    @test isapprox(response, 0.056716; atol=1e-3)

    test_pupadurationrossi = PupaDurationRossi(0.027, -1.7, 27.7)
    response = get_temperature_response(27.0, test_pupadurationrossi)
    @test isapprox(response, 1.483000; atol=1e-3)

    test_adultmortalityrossi = AdultMortalityRossi(0.053, 0.081, 23.0, 6.375 * 10^(-4))
    response = get_temperature_response(27.0, test_adultmortalityrossi)
    @test isapprox(response, 0.024702; atol=1e-3)

    # AedesAegypti: El Moustaid et al 2019
    test_eggdurationmoustaid = EggDurationMoustaid(87.1722, 1.1575, 0.0136)
    response = get_temperature_response(27.0, test_eggdurationmoustaid)
    @test isapprox(response, 0.015189; atol=1e-3)

    test_eggmortalitymoustaid =
        EggMortalityMoustaid(30.33, 9.395, 3.6576, 110.0014, 0.00077)
    response = get_temperature_response(27.0, test_eggmortalitymoustaid, 0.015189)
    @test isapprox(response, 0.005044; atol=1e-3)

    test_larvadurationmoustaid = LarvaDurationMoustaid(29.0889, 1.9571, 0.03357)
    response = get_temperature_response(27.0, test_larvadurationmoustaid)
    @test isapprox(response, 1.389409; atol=1e-3)

    test_larvamortalitymoustaid =
        LarvaMortalityMoustaid(36.79, 2.354, 2.01938, 73.1928, 0.00046)
    response = get_temperature_response(27.0, test_larvamortalitymoustaid, 1.389409)
    @test isapprox(response, 0.058680; atol=1e-3)

    test_pupadurationmoustaid = PupaDurationMoustaid(20.5074, 1.0954, 0.0153)
    response = get_temperature_response(27.0, test_pupadurationmoustaid)
    @test isapprox(response, 0.479547; atol=1e-3)

    test_pupamortalitymoustaid =
        PupaMortalityMoustaid(37.4679, 9.96278, 1.0, 38.0, 0.03487954)
    response = get_temperature_response(27.0, test_pupamortalitymoustaid, 0.479547)
    @test isapprox(response, 0.0; atol=1e-3)

    test_adultmortalitymoustaid = AdultMortalityMoustaid(37.73, 9.16, -0.149)
    response = get_temperature_response(27.0, test_adultmortalitymoustaid)
    @test isapprox(response, 0.035060; atol=1e-3)

    # AnophelesGambiae: Abiodun et al 2016
    test_eggdurationabiodun =
        EggDurationAbiodun(6.0, 1.011, 20.212, 1.0, 12.096, 4.839, 1.0)
    response = get_temperature_response(27.0, test_eggdurationabiodun)
    @test isapprox(response, 1.166976; atol=1e-3)

    test_eggmortalityabiodun =
        EggMortalityAbiodun(6.0, 1.011, 20.212, 1.0, 12.096, 4.839, 1.0)
    response = get_temperature_response(27.0, test_eggmortalityabiodun)
    @test isapprox(response, 0.424469; atol=1e-3)

    test_larvadurationabiodun =
        LarvaDurationAbiodun(6.0, 8.130, 13.794, 1.0, 12.096, 4.839, 1.0, 1.011, 20.212)
    response = get_temperature_response(27.0, test_larvadurationabiodun)
    @test isapprox(response, 7.069472; atol=1e-3)

    test_larvamortalityabiodun =
        LarvaMortalityAbiodun(6.0, 8.130, 13.794, 1.0, 12.096, 4.839, 1.0, 1.011, 20.212)
    response = get_temperature_response(27.0, test_larvamortalityabiodun)
    @test isapprox(response, 0.885669; atol=1e-3)

    test_pupadurationabiodun = PupaDurationAbiodun(
        6.0,
        8.560,
        20.654,
        1.0,
        19.759,
        6.827,
        1.0,
        8.130,
        13.794,
        12.096,
        4.839,
    )
    response = get_temperature_response(27.0, test_pupadurationabiodun, 8.236448)
    @test isapprox(response, 0.928051; atol=1e-3)

    test_pupamortalityabiodun = PupaMortalityAbiodun(
        6.0,
        8.560,
        20.654,
        1.0,
        19.759,
        6.827,
        1.0,
        8.130,
        13.794,
        12.096,
        4.839,
    )
    response = get_temperature_response(27.0, test_pupamortalityabiodun)
    @test isapprox(response, 0.896625; atol=1e-3)

    test_adultmortalityabiodun = AdultMortalityAbiodun(4.4, 1.31, 0.03)
    response = get_temperature_response(27.0, test_adultmortalityabiodun)
    @test isapprox(response, 0.109890; atol=1e-3)
end
