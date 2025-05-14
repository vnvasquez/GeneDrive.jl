@testset "Decision Model" begin

    # Shared
    i = JuMP.optimizer_with_attributes(Ipopt.Optimizer)
    GAE = JuMP.GenericAffExpr{Float64, VariableRef}
    GQE = JuMP.GenericQuadExpr{Float64, VariableRef}

    ## Constraints: M
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
    test_strategy = [1 => test_constraints]
    test_node_strategy = NodeStrategy(1, test_strategy)

    test_prob = create_decision_model(
        test_node,
        test_short_tspan;
        node_strategy=test_node_strategy,
        node_species=[test_species], # now requires vector
        optimizer=i,
    )

    @test JuMP.num_variables(test_prob) == 460
    @test JuMP.num_constraints(test_prob, GAE, MOI.EqualTo{Float64}) == 150
    @test JuMP.num_constraints(test_prob, GAE, MOI.GreaterThan{Float64}) == 30
    @test JuMP.num_constraints(test_prob, GAE, MOI.LessThan{Float64}) == 61
    @test JuMP.num_constraints(test_prob, GQE, MOI.EqualTo{Float64}) == 90
    @test JuMP.num_constraints(test_prob, VariableRef, MOI.EqualTo{Float64}) == 90
    @test JuMP.num_constraints(test_prob, VariableRef, MOI.GreaterThan{Float64}) == 370
    @test JuMP.num_constraints(test_prob, VariableRef, MOI.LessThan{Float64}) == 10

    ## Constraints: F && slack_small && _create_default_solver
    test_constraints = ReleaseStrategy(
        release_this_gene_index=test_release_gene,
        release_this_life_stage=Female,
    )
    test_strategy = [1 => test_constraints]
    test_node_strategy = NodeStrategy(1, test_strategy)

    test_prob = create_decision_model(
        test_node,
        test_short_tspan;
        node_strategy=test_node_strategy,
        node_species=[test_species],    # now requires vector
        slack_small=true,
    )

    @test JuMP.num_variables(test_prob) == 463
    @test JuMP.num_constraints(test_prob, GAE, MOI.EqualTo{Float64}) == 150
    @test JuMP.num_constraints(test_prob, GAE, MOI.GreaterThan{Float64}) == 90
    @test JuMP.num_constraints(test_prob, GAE, MOI.LessThan{Float64}) == 121
    @test JuMP.num_constraints(test_prob, GQE, MOI.EqualTo{Float64}) == 90
    @test JuMP.num_constraints(test_prob, VariableRef, MOI.EqualTo{Float64}) == 80
    @test JuMP.num_constraints(test_prob, VariableRef, MOI.GreaterThan{Float64}) == 383
    @test JuMP.num_constraints(test_prob, VariableRef, MOI.LessThan{Float64}) == 13

    ## Constraints: Mixed && slack_large
    test_constraints = ReleaseStrategy(
        release_this_gene_index=test_release_gene,
        release_this_life_stage=(Male, Female),
    )
    test_strategy = [1 => test_constraints]
    test_node_strategy = NodeStrategy(1, test_strategy)

    test_prob = create_decision_model(
        test_node,
        test_short_tspan;
        node_strategy=test_node_strategy,
        node_species=[test_species],
        optimizer=i,
        slack_large=true,
    )

    @test JuMP.num_variables(test_prob) == 490
    @test JuMP.num_constraints(test_prob, GAE, MOI.EqualTo{Float64}) == 160
    @test JuMP.num_constraints(test_prob, GAE, MOI.GreaterThan{Float64}) == 120
    @test JuMP.num_constraints(test_prob, GAE, MOI.LessThan{Float64}) == 152
    @test JuMP.num_constraints(test_prob, GQE, MOI.EqualTo{Float64}) == 90
    @test JuMP.num_constraints(test_prob, VariableRef, MOI.EqualTo{Float64}) == 60
    @test JuMP.num_constraints(test_prob, VariableRef, MOI.GreaterThan{Float64}) == 430
    @test JuMP.num_constraints(test_prob, VariableRef, MOI.LessThan{Float64}) == 40
end
