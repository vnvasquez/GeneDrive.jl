@testset "Decision Model" begin
    test_species = AnophelesGambiae
    test_anoph = update_population_size(stages_abiodun(), 500)
    test_mendelian = make_organisms(test_species, genetics_mendelian(), test_anoph)
    test_short_timeseries = TimeSeriesTemperature([
        27.0,
        27.5,
        28.05,
        28.05,
        28.3,
        27.8,
        27.9,
        28.5,
        27.75,
        26.8,
    ])
    test_short_tspan = (1, length(test_short_timeseries.values))
    test_node = Node(:TestNode, test_mendelian, test_short_timeseries, (1.0, 1.0))

    test_wild_gene = get_wildtype(test_node, test_species)
    test_release_gene = get_homozygous_modified(test_node, test_species)
    test_constraints = ReleaseStrategy(
        release_this_gene_index=test_release_gene,
        release_this_life_stage=Male,
    )
    test_strategy = Dict(1 => test_constraints)

    test_prob = create_decision_model(
        test_node,
        test_short_tspan;
        node_strategy=test_strategy,
        species=test_species,
        optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer),
        slack_small=true,
    )

    GAE = JuMP.GenericAffExpr{Float64, VariableRef}
    GQE = JuMP.GenericQuadExpr{Float64, VariableRef}
    @test JuMP.num_variables(test_prob) == 463
    @test JuMP.num_constraints(test_prob, GAE, MOI.EqualTo{Float64}) == 150
    @test JuMP.num_constraints(test_prob, GAE, MOI.GreaterThan{Float64}) == 30
    @test JuMP.num_constraints(test_prob, GAE, MOI.LessThan{Float64}) == 61
    @test JuMP.num_constraints(test_prob, GQE, MOI.EqualTo{Float64}) == 90
    @test JuMP.num_constraints(test_prob, VariableRef, MOI.EqualTo{Float64}) == 90
    @test JuMP.num_constraints(test_prob, VariableRef, MOI.GreaterThan{Float64}) == 373
    @test JuMP.num_constraints(test_prob, VariableRef, MOI.LessThan{Float64}) == 13
end
