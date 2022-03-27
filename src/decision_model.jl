################################################################################
#                                  Create Model                                 #
################################################################################
"""
    create_decision_model(network::Network, tspan; node_strategy::Dict, do_binary::Bool=false)

Build mathematical program. Problem created as an NLP (do_binary=false) or MINLP (do_binary=true).
"""
function create_decision_model(
    network::Network,
    tspan;
    node_strategy::Dict,
    do_binary::Bool=false,
)
    for (key_node, node) in enumerate(values(get_nodes(network)))
        if length(collect(tspan[1]:tspan[2])) !== length(node.temperature.values)
            @warn(
                "Temperature timeseries in node $key_node does not match `tspan` length. Problem will fail or result will be incorrect."
            )
        end
    end

    ##################
    # Initialization
    ##################
    initial_condition, density_net = init_network!(network)

    ##################
    # Solver(s)
    ##################
    ipopt_def = JuMP.optimizer_with_attributes(
        Ipopt.Optimizer,
        "linear_solver" => "pardiso",
        "print_level" => 1,
    )

    i = JuMP.optimizer_with_attributes(
        Juniper.Optimizer,
        "nl_solver" => ipopt_def,
        "mip_solver" => JuMP.optimizer_with_attributes(Gurobi.Optimizer),
    )

    ##################
    #  Model Creation
    ##################
    if do_binary
        model = JuMP.Model(i)
        @info(@info("Ensure that the solver(s) being called are installed: $(i)"))
    else
        model = JuMP.Model(ipopt_def)
        @info("Ensure that the solver(s) being called are installed: $(ipopt_def)")
    end

    ##################
    # Parameters
    ##################
    data = _optimization_info(network, tspan)

    # TODO: fix
    nE = data[1]["organism"][1]["substage_count"][Egg]
    nL = data[1]["organism"][1]["substage_count"][Larva]
    nP = data[1]["organism"][1]["substage_count"][Pupa]
    nM = data[1]["organism"][1]["substage_count"][Male]
    nF = data[1]["organism"][1]["substage_count"][Female]
    gene_count = data[1]["organism"][1]["gene_count"]
    homozygous_modified = data[1]["organism"][1]["homozygous_modified"]
    wildtype = data[1]["organism"][1]["wildtype"]
    species = AedesAegypti

    # TODO: fix
    densE = data[1]["organism"][1]["stage_density"][Egg]
    densL = data[1]["organism"][1]["stage_density"][Larva]
    densP = data[1]["organism"][1]["stage_density"][Pupa]
    densM = data[1]["organism"][1]["stage_density"][Male]
    densF = data[1]["organism"][1]["stage_density"][Female]

    # SETS
    ##################
    T = 1:tspan[end]
    N = 1:data["total_node_count"]
    O = 1:data[1]["node_organism_count"] #TODO: fix

    # Stage/substage sets TODO: fix
    SE = 1:nE
    SL = 1:nL
    SP = 1:nP
    SM = 1:nM
    SF = 1:nF

    # Initial conditions mapped to stage/substage sets TODO: fix
    SE_map = SE
    SL_map = (nE + 1):(nE + nL)
    SP_map = (nE + nL + 1):(nE + nL + nP)
    SM_map = nE + nL + nP + 1
    SF_map = (nE + nL + nP + nM + 1):(nE + nL + nP + nM + nF)

    # Genes TODO: fix
    G = 1:gene_count

    # Add sets to model object TODO: fix
    model.obj_dict[:Sets] = Dict(
        :N => N,
        :O => O,
        :SE => SE,
        :SL => SL,
        :SP => SP,
        :SM => SM,
        :SF => SF,
        :G => G,
        :T => T,
    )

    # DECLARE VARIABLES_1: Life Stages
    ###########################################
    JuMP.@variable(model, E[N, O, SE, G, T] >= 0)
    JuMP.@variable(model, L[N, O, SL, G, T] >= 0)
    JuMP.@variable(model, P[N, O, SP, G, T] >= 0)
    JuMP.@variable(model, M[N, O, SM, G, T] >= 0)
    JuMP.@variable(model, F[N, O, SF, G, T] >= 0)

    # DECLARE VARIABLES_2: Binary nodes
    ###########################################
    if do_binary
        JuMP.@variable(model, release_location[N], Bin)
    else
        JuMP.@variable(model, release_location[N])
        JuMP.fix.(release_location, 1.0)
    end

    # DECLARE VARIABLES_3: Controls
    ###########################################
    JuMP.@variable(model, 0.0 <= control_M[N, O, SM, G, T])
    JuMP.@variable(model, 0.0 <= control_F[N, O, SF, G, T])

    # WARMSTART VARIABLES_1: Lifestages
    ###########################################
    stages = [Egg, Larva, Pupa, Male, Female]
    initialcond_dict = Dict(
        s => Dict(n => Dict{Int, Any}(o => nothing for o in O) for n in N) for s in stages
    )
    for (ix, node_name) in enumerate(N)
        for (jx, organism) in enumerate(O)
            initialcond_dict[Egg][node_name][organism] =
                initial_condition.x[ix].x[jx][SE, :]
            initialcond_dict[Larva][node_name][organism] =
                initial_condition.x[ix].x[jx][(nE + 1):(nL + nE), :]
            initialcond_dict[Pupa][node_name][organism] =
                initial_condition.x[ix].x[jx][(nE + nL + 1):(nE + nL + nP), :]
            initialcond_dict[Male][node_name][organism] =
                initial_condition.x[ix].x[jx][nE + nL + nP + 1, :]'
            initialcond_dict[Female][node_name][organism] =
                initial_condition.x[ix].x[jx][(nE + nL + nP + 2):end, :]
        end
    end
    for node_name in N, organism in O, t in T
        JuMP.set_start_value.(
            E[node_name, organism, :, :, t].data,
            initialcond_dict[Egg][node_name][organism],
        )
        JuMP.set_start_value.(
            L[node_name, organism, :, :, t].data,
            initialcond_dict[Larva][node_name][organism],
        )
        JuMP.set_start_value.(
            P[node_name, organism, :, :, t].data,
            initialcond_dict[Pupa][node_name][organism],
        )
        JuMP.set_start_value.(
            M[node_name, organism, :, :, t].data,
            initialcond_dict[Male][node_name][organism],
        )
        JuMP.set_start_value.(
            F[node_name, organism, :, :, t].data,
            initialcond_dict[Female][node_name][organism],
        )
    end

    # EXPRESSIONS: Migration
    ###########################################
    A = get_migration(network, species)
    JuMP.@expression(
        model,
        migration_E[n in N, o in O, s in SE, g in G, t in T],
        A[SE_map[s], g][n, :]' * E[:, o, s, g, t]
    )
    JuMP.@expression(
        model,
        migration_L[n in N, o in O, s in SL, g in G, t in T],
        A[SL_map[s], g][n, :]' * L[:, o, s, g, t]
    )
    JuMP.@expression(
        model,
        migration_P[n in N, o in O, s in SP, g in G, t in T],
        A[SP_map[s], g][n, :]' * P[:, o, s, g, t]
    )
    JuMP.@expression(
        model,
        migration_M[n in N, o in O, s in SM, g in G, t in T],
        A[SM_map[s], g][n, :]' * M[:, o, s, g, t]
    )

    # CONSTRAINTS_A: Life Stages
    ###########################################

    #### EGGS
    JuMP.@constraint(
        model,
        E_con_A0[n in N, o in O, s in [SE[1]], g in G, t in [T[1]]],
        E[n, o, s, g, t] ==
        initialcond_dict[Egg][n][o][SE_map[s], g] +
        sum(
            (data[n]["organism"][o]["genetics"].likelihood[
                :,
                :,
                g,
            ] .* data[n]["organism"][o]["genetics"].Τ[
                :,
                :,
                g,
            ] .* data[n]["organism"][o]["genetics"].S[g] * data[n]["organism"][o]["genetics"].Β[g] .* F[
                n,
                o,
                :,
                :,
                t,
            ].data)[
                :,
                :,
            ],
        ) -
        E[n, o, s, g, t] * (
            data[n]["organism"][o]["stage_temperature_response"][Egg][n][o][g][t][1] *
            compute_density(densE, sum(E[n, o, :, :, t])) +
            data[n]["organism"][o]["stage_temperature_response"][Egg][n][o][g][t][2] * nE
        ) + migration_E[n, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        E_con_A1[n in N, o in O, s in [SE[1]], g in G, t in T[2:end]],
        E[n, o, s, g, t] ==
        E[n, o, s, g, t - 1] +
        sum(
            (data[n]["organism"][o]["genetics"].likelihood[
                :,
                :,
                g,
            ] .* data[n]["organism"][o]["genetics"].Τ[
                :,
                :,
                g,
            ] .* data[n]["organism"][o]["genetics"].S[g] * data[n]["organism"][o]["genetics"].Β[g] .* F[
                n,
                o,
                :,
                :,
                t,
            ].data)[
                :,
                :,
            ],
        ) -
        E[n, o, s, g, t] * (
            data[n]["organism"][o]["stage_temperature_response"][Egg][n][o][g][t][1] *
            compute_density(densE, sum(E[n, o, :, :, t])) +
            data[n]["organism"][o]["stage_temperature_response"][Egg][n][o][g][t][2] * nE
        ) + migration_E[n, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        E_con_B0[n in N, o in O, s in SE[2:end], g in G, t in [T[1]]],
        E[n, o, s, g, t] ==
        initial_condition.x[n].x[o][SE_map[s], g] +
        data[n]["organism"][o]["stage_temperature_response"][Egg][n][o][g][t][2] *
        nE *
        E[n, o, s - 1, g, t] -
        E[n, o, s, g, t] * (
            data[n]["organism"][o]["stage_temperature_response"][Egg][n][o][g][t][1] *
            compute_density(densE, sum(E[n, o, :, :, t])) +
            data[n]["organism"][o]["stage_temperature_response"][Egg][n][o][g][t][2] * nE
        ) + migration_E[n, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        E_con_B1[n in N, o in O, s in SE[2:end], g in G, t in T[2:end]],
        E[n, o, s, g, t] ==
        E[n, o, s, g, t - 1] +
        data[n]["organism"][o]["stage_temperature_response"][Egg][n][o][g][t][2] *
        nE *
        E[n, o, s - 1, g, t] -
        E[n, o, s, g, t] * (
            data[n]["organism"][o]["stage_temperature_response"][Egg][n][o][g][t][1] *
            compute_density(densE, sum(E[n, o, :, :, t])) +
            data[n]["organism"][o]["stage_temperature_response"][Egg][n][o][g][t][2] * nE
        ) + migration_E[n, o, s, g, t]
    )

    #### LARVAE
    JuMP.@constraint(
        model,
        L_con_A0[n in N, o in O, s in [SL[1]], g in G, t in [T[1]]],
        L[n, o, s, g, t] ==
        initial_condition.x[n].x[o][SL_map[s], g] +
        data[n]["organism"][o]["stage_temperature_response"][Egg][n][o][g][t][2] *
        nE *
        E[n, o, end, g, t] -
        L[n, o, s, g, t] * (
            data[n]["organism"][o]["stage_temperature_response"][Larva][n][o][g][t][1] *
            compute_density(densL, sum(L[n, o, :, :, t])) +
            data[n]["organism"][o]["stage_temperature_response"][Larva][n][o][g][t][2] * nL
        ) + migration_L[n, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        L_con_A1[n in N, o in O, s in [SL[1]], g in G, t in T[2:end]],
        L[n, o, s, g, t] ==
        L[n, o, s, g, t - 1] +
        data[n]["organism"][o]["stage_temperature_response"][Egg][n][o][g][t][2] *
        nE *
        E[n, o, end, g, t] -
        L[n, o, s, g, t] * (
            data[n]["organism"][o]["stage_temperature_response"][Larva][n][o][g][t][1] *
            compute_density(densL, sum(L[n, o, :, :, t])) +
            data[n]["organism"][o]["stage_temperature_response"][Larva][n][o][g][t][2] * nL
        ) + migration_L[n, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        L_con_B0[n in N, o in O, s in SL[2:end], g in G, t in [T[1]]],
        L[n, o, s, g, t] ==
        initial_condition.x[n].x[o][SL_map[s], g] +
        data[n]["organism"][o]["stage_temperature_response"][Larva][n][o][g][t][2] *
        nL *
        L[n, o, s - 1, g, t] -
        L[n, o, s, g, t] * (
            data[n]["organism"][o]["stage_temperature_response"][Larva][n][o][g][t][1] *
            compute_density(densL, sum(L[n, o, :, :, t])) +
            data[n]["organism"][o]["stage_temperature_response"][Larva][n][o][g][t][2] * nL
        ) + migration_L[n, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        L_con_B1[n in N, o in O, s in SL[2:end], g in G, t in T[2:end]],
        L[n, o, s, g, t] ==
        L[n, o, s, g, t - 1] +
        data[n]["organism"][o]["stage_temperature_response"][Larva][n][o][g][t][2] *
        nL *
        L[n, o, s - 1, g, t] -
        L[n, o, s, g, t] * (
            data[n]["organism"][o]["stage_temperature_response"][Larva][n][o][g][t][1] *
            compute_density(densL, sum(L[n, o, :, :, t])) +
            data[n]["organism"][o]["stage_temperature_response"][Larva][n][o][g][t][2] * nL
        ) + migration_L[n, o, s, g, t]
    )

    #### PUPAE
    JuMP.@constraint(
        model,
        P_con_A0[n in N, o in O, s in [SP[1]], g in G, t in [T[1]]],
        P[n, o, s, g, t] ==
        initial_condition.x[n].x[o][SP_map[s], g] +
        data[n]["organism"][o]["stage_temperature_response"][Larva][n][o][g][t][2] *
        nL *
        L[n, o, end, g, t] -
        P[n, o, s, g, t] * (
            data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][1] *
            compute_density(densP, sum(P[n, o, :, :, t])) +
            data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][2] * nP
        ) + migration_P[n, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        P_con_A1[n in N, o in O, s in [SP[1]], g in G, t in T[2:end]],
        P[n, o, s, g, t] ==
        P[n, o, s, g, t - 1] +
        data[n]["organism"][o]["stage_temperature_response"][Larva][n][o][g][t][2] *
        nL *
        L[n, o, end, g, t] -
        P[n, o, s, g, t] * (
            data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][1] *
            compute_density(densP, sum(P[n, o, :, :, t])) +
            data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][2] * nP
        ) + migration_P[n, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        P_con_B0[n in N, o in O, s in SP[2:end], g in G, t in [T[1]]],
        P[n, o, s, g, t] ==
        initial_condition.x[n].x[o][SP_map[s], g] +
        data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][2] *
        nP *
        P[n, o, s - 1, g, t] -
        P[n, o, s, g, t] * (
            data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][1] *
            compute_density(densP, sum(P[n, o, :, :, t])) +
            data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][2] * nP
        ) + migration_P[n, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        P_con_B1[n in N, o in O, s in SP[2:end], g in G, t in T[2:end]],
        P[n, o, s, g, t] ==
        P[n, o, s, g, t - 1] +
        data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][2] *
        nP *
        P[n, o, s - 1, g, t] -
        P[n, o, s, g, t] * (
            data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][1] *
            compute_density(densP, sum(P[n, o, :, :, t])) +
            data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][2] * nP
        ) + migration_P[n, o, s, g, t]
    )

    #### MALES
    JuMP.@constraint(
        model,
        M_con_0[n in N, o in O, s in SM, g in G, t in [T[1]]],
        M[n, o, s, g, t] ==
        initial_condition.x[n].x[o][SM_map[s], g] +
        (1 - data[n]["organism"][o]["genetics"].Φ[g]) *
        data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][2] *
        nP *
        P[n, o, end, g, t] *
        data[n]["organism"][o]["genetics"].Ξ_m[g] -
        (1 + data[n]["organism"][o]["genetics"].Ω[g]) *
        data[n]["organism"][o]["stage_temperature_response"][Male][n][o][g][t][1] *
        M[n, o, s, g, t] *
        compute_density(densM, sum(M[n, o, :, :, t])) + migration_M[n, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        M_con_1[n in N, o in O, s in SM, g in G, t in T[2:end]],
        M[n, o, s, g, t] ==
        M[n, o, s, g, t - 1] +
        (1 - data[n]["organism"][o]["genetics"].Φ[g]) *
        data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][2] *
        nP *
        P[n, o, end, g, t] *
        data[n]["organism"][o]["genetics"].Ξ_m[g] -
        (1 + data[n]["organism"][o]["genetics"].Ω[g]) *
        data[n]["organism"][o]["stage_temperature_response"][Male][n][o][g][t][1] *
        M[n, o, s, g, t] *
        compute_density(densM, sum(M[n, o, :, :, t])) +
        control_M[n, o, s, g, t] +
        migration_M[n, o, s, g, t]
    )

    #### MATING
    JuMP.@constraint(
        model,
        mate_bound[n in N, o in O, g in G, t in T],
        M[n, o, 1, g, t] * data[n]["organism"][o]["genetics"].Η[g] <=
        (sum(M[n, o, 1, i, t] * data[n]["organism"][o]["genetics"].Η[i] for i in G))
    )

    #### FEMALES:
    JuMP.@NLconstraint(
        model,
        F_con_0[n in N, o in O, s in SF, g in G, t in [T[1]]],
        F[n, o, s, g, t] ==
        initial_condition.x[n].x[o][SF_map[s], g] +
        (
            M[n, o, 1, g, t] * data[n]["organism"][o]["genetics"].Η[g] / (
                1e-6 +
                sum(M[n, o, 1, i, t] * data[n]["organism"][o]["genetics"].Η[i] for i in G)
            )
        ) * (
            data[n]["organism"][o]["genetics"].Φ[s] *
            data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][2] *
            nP *
            P[n, o, end, s, t] *
            data[n]["organism"][o]["genetics"].Ξ_f[s]
        ) -
        (1 + data[n]["organism"][o]["genetics"].Ω[g]) *
        data[n]["organism"][o]["stage_temperature_response"][Female][n][o][g][t][1] *
        F[n, o, s, g, t] + sum(A[SF_map[s], g][n, i] * F[i, o, s, g, t] for i in N)
    )

    JuMP.@NLconstraint(
        model,
        F_con_1[n in N, o in O, s in SF, g in G, t in T[2:end]],
        F[n, o, s, g, t] ==
        F[n, o, s, g, t - 1] +
        (
            M[n, o, 1, g, t] * data[n]["organism"][o]["genetics"].Η[g] / (
                1e-6 +
                sum(M[n, o, 1, i, t] * data[n]["organism"][o]["genetics"].Η[i] for i in G)
            )
        ) * (
            data[n]["organism"][o]["genetics"].Φ[s] *
            data[n]["organism"][o]["stage_temperature_response"][Pupa][n][o][g][t][2] *
            nP *
            P[n, o, end, s, t] *
            data[n]["organism"][o]["genetics"].Ξ_f[s]
        ) -
        (1 + data[n]["organism"][o]["genetics"].Ω[g]) *
        data[n]["organism"][o]["stage_temperature_response"][Female][n][o][g][t][1] *
        F[n, o, s, g, t] +
        control_F[n, o, s, g, t] +
        sum(A[SF_map[s], g][n, i] * F[i, o, s, g, t] for i in N)
    )

    # CONSTRAINTS_B: Controls
    ###########################################
    _calculate_release_constraints(
        network,
        tspan,
        homozygous_modified,
        wildtype,
        do_binary,
        model,
        data,
        node_strategy,
    )

    return model
end

################################################################################
#                                  Solve Model                                 #
################################################################################
"""
    solve_decision_model(model::JuMP.Model, objective_function::Nothing=nothing; kwargs...)

Solve mathematical program using default objective function `@objective(model, Min, 0)`. This permits comparison to dynamic population model output.
"""
function solve_decision_model(
    model::JuMP.Model,
    objective_function::Nothing=nothing;
    kwargs...,
)

    #Re-set constraints by fixing controls to zero
    control_M = model[:control_M]
    JuMP.fix.(control_M, 0.0; force=true)
    control_F = model[:control_F]
    JuMP.fix.(control_F, 0.0; force=true)

    # Permit optimizer to act as a nonlinear solver
    JuMP.@objective(model, Min, 0)
    JuMP.optimize!(model)
    @info("No objective function specified. Default supplied: `@objective(model, Min, 0)`")

    return model
end
