@testset "Dynamic Model" begin

    # Shared
    test_solver = OrdinaryDiffEq.Tsit5()
    test_tspan = (1, 365)

    # Temperature: Timeseries, Genetics: Mendelian, Species: AnophelesGambiae, Releases: no, Shocks: no, Network: no
    test_anoph = update_population_size(stages_abiodun(), 500)
    test_mendelian = make_organisms(AnophelesGambiae, genetics_mendelian(), test_anoph)
    test_node = Node(:TestNode, test_mendelian, example_temperature_timeseries, (1.0, 1.0))
    test_sol = solve_dynamic_model(test_node, test_solver, test_tspan)
    plot_dynamic_mendelian_females(test_node, test_sol)
    @test test_sol.retcode == :Success
    # @test LinearAlgebra.norm(realanswer - testanswer) < 1e-3

    ### Temperature: Sinusoid, Genetics: RIDL, Species: Aedes, Releases: yes, Shocks: yes, Network: no
    test_aedes = update_population_size(stages_rossi(), 500)
    test_ridl = make_organisms(AedesAegypti, genetics_ridl(), test_aedes)
    test_node = Node(:TestNode, test_ridl, example_temperature_sinusoidal, (1.0, 1.0))
    test_wild_gene = get_wildtype(test_node, AedesAegypti)
    test_release_gene = get_homozygous_modified(test_node, AedesAegypti)
    test_release_size = 100
    test_release_times = [5.0, 10.0]
    test_release = Release(
        test_node,
        AedesAegypti,
        Male,
        test_release_gene,
        test_release_times,
        test_release_size,
    )
    test_shocks = TemperatureShockData(test_node, [(20.0, 30.0)], [2.0])
    test_sol =
        solve_dynamic_model(test_node, [test_release], test_shocks, test_solver, test_tspan)
    plot_dynamic_ridl_females(test_node, test_sol)
    @test test_sol.retcode == :Success

    ## Temperature: Constant, Genetics: HGD, Species: Aedes, Releases: no, Shocks: no, Network: no
    test_aedes = update_population_size(stages_rossi(), 500)
    test_hgd = make_organisms(AedesAegypti, genetics_mcr(), test_aedes)
    test_node = Node(:TestNode, test_hgd, example_temperature_constant, (1.0, 1.0))
    test_sol = solve_dynamic_model(test_node, test_solver, test_tspan)
    plot_dynamic_mcr_females(test_node, test_sol)
    @test test_sol.retcode == :Success

    ## Temperature: Constant, Genetics: Wolbachia, Species: Aedes, Releases: no, Shocks: no, Network: yes
    test_aedes = update_population_size(stages_rossi(), 500)
    test_wolb = make_organisms(AedesAegypti, genetics_wolbachia(), test_aedes)
    test_node = Node(:TestNode, test_wolb, example_temperature_constant, (1.0, 1.0))
    test_node2 = Node(:TestNode2, test_wolb, example_temperature_constant, (2.0, 2.0))
    test_net = Network(:TestNet, test_node, test_node2)
    move_rate = 0.002
    test_migration = Dict(
        ("Male", "WW") => Dict(
            (:TestNode, :TestNode2) => move_rate,
            (:TestNode2, :TestNode) => move_rate,
        ),
        ("Male", "ww") => Dict(
            (:TestNode, :TestNode2) => move_rate,
            (:TestNode2, :TestNode) => move_rate,
        ),
        ("Female", "WW") => Dict(
            (:TestNode, :TestNode2) => move_rate,
            (:TestNode2, :TestNode) => move_rate,
        ),
        ("Female", "ww") => Dict(
            (:TestNode, :TestNode2) => move_rate,
            (:TestNode2, :TestNode) => move_rate,
        ),
    )
    assign_migration!(test_net, test_migration, AedesAegypti)
    test_sol = solve_dynamic_model(test_net, test_solver, test_tspan)
    @test test_sol.retcode == :Success
end
