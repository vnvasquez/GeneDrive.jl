################################################################################
#                                  Create Model                                 #
################################################################################

function _create_default_solvers(do_binary::Bool)
    ipopt_def = JuMP.optimizer_with_attributes(Ipopt.Optimizer)

    if !do_binary
        return ipopt_def
    end

    i = JuMP.optimizer_with_attributes(
        Juniper.Optimizer,
        "nl_solver" => ipopt_def,
        "mip_solver" => JuMP.optimizer_with_attributes(Cbc.Optimizer),
    )

    return i
end

"""
    create_decision_model(network::Network, tspan; node_strategy::DataStructures.OrderedDict, species::Type{<:Species}=AedesAegypti,do_binary::Bool=false, optimizer=nothing)

Build mathematical program. Problem created as an NLP (do_binary=false) or MINLP (do_binary=true).

# Arguments

  - `network`: Data for one or more interconnected spatial nodes.
  - `tspan`: Time horizon of the problem, defined in days.
  - `node_strategy`: Operational constraints specific to organism types within each node. 
  - `node_species`: Species present in the node, defined as a vector. 
  - `do_binary`: Boolean indicating whether the problem is formulated as a MINLP (true) or NLP (false); defaults to false.
  - `optimizer`: Optimization solver to be used.
  - `slack_small`: Boolean flag to enable small slack variables within selected constraints. These introduce minor constraint relaxations with corresponding penalties in the objective function. Defaults to false.
  - `slack_large`: Boolean flag to enable large slack variables for constraints where greater flexibility may be required to preserve model feasibility. Defaults to false.
"""
function create_decision_model(
    network::Network,
    tspan;
    node_strategy,
    node_species, 
    do_binary::Bool=false,
    optimizer=nothing,
    slack_small=false,
    slack_large=false,
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
    # Solver
    ##################
    i = optimizer === nothing ? _create_default_solvers(do_binary) : optimizer

    ##################
    # Model
    ##################
    model = JuMP.Model(i)

    ##################
    # Parameters
    ##################
    data = _optimization_info(network, tspan)
    organism_data = data[1]["scenario"][1]
    gene_count = organism_data["organism"][1]["gene_count"]

    # TODO: fix
    nE = organism_data["organism"][1]["substage_count"][Egg]
    nL = organism_data["organism"][1]["substage_count"][Larva]
    nP = organism_data["organism"][1]["substage_count"][Pupa]
    nM = organism_data["organism"][1]["substage_count"][Male]
    nF = organism_data["organism"][1]["substage_count"][Female]

    # TODO: fix
    densE = organism_data["organism"][1]["stage_density"][Egg]
    densL = organism_data["organism"][1]["stage_density"][Larva]
    densP = organism_data["organism"][1]["stage_density"][Pupa]
    densM = organism_data["organism"][1]["stage_density"][Male]
    densF = organism_data["organism"][1]["stage_density"][Female]

    # SETS
    ##################
    T = 1:tspan[end]
    N = 1:data["total_node_count"]
    C = 1:data["total_scenario_count"]
    O = 1:organism_data["node_organism_count"] 
    G = 1:gene_count

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
 
    #= Genes
    G = Dict()
    for org_key in keys(organism_data["organism"])
        gene_count = organism_data["organism"][org_key]["gene_count"]
        G[org_key] = 1:gene_count
    end 
    =#

    # Add sets to model object TODO: fix
    model.obj_dict[:Sets] = Dict(
        :N => N,
        :C => C,
        :O => O,
        :SE => SE,
        :SL => SL,
        :SP => SP,
        :SM => SM,
        :SF => SF,
        :G => G,
        :T => T,
    )

    # Add probabilities to model ext (because not optimizing, using as extension)
    model.ext[:Probabilities] =
        Dict(n => Dict{Int, Float64}(c => 0.0 for c in C) for n in N)

    #for (cx, prob) in enumerate(get_probability(node.temperature))
    for (node_number, node) in model.ext[:Probabilities]
        for cx in keys(node)
            node[cx] = data[node_number]["scenario"][cx]["probability"]
        end
    end
  #=
    for n in N, o in O, s in SE, t in T
        for g in G
            JuMP.set_start_value(E[n, o ,s,g, t], ini_cond[....])
        end
    end
    =# 
    # G = 1:3
    # DECLARE VARIABLES_1: Life Stages
    ###########################################
    #JuMP.@variable(model, E[N, C, o in O, SE, g in G[o], T] >= 0) 
    JuMP.@variable(model, E[N, C, O, SE, G, T] >= 0) # maybe?
    #error()
    JuMP.@variable(model, L[N, C, O, SL, G, T] >= 0)
    JuMP.@variable(model, P[N, C, O, SP, G, T] >= 0)
    JuMP.@variable(model, M[N, C, O, SM, G, T] >= 0)
    JuMP.@variable(model, F[N, C, O, SF, G, T] >= 0)
    #JuMP.@variable(model, P_slack_negative[N, O, G, T] == 0)
    if slack_small
        JuMP.@variable(model, 0 <= P_slack_positive[N, C, O, G, [1]] <= 1)
    elseif slack_large
        JuMP.@variable(model, 0 <= P_slack_positive[N, C, O, G, T] <= 1)
    else
        @info("No slack specified.")
    end

    # DECLARE VARIABLES_2: Controls
    ###########################################
    JuMP.@variable(model, 0.0 <= control_M[N, O, SM, G, T])
    JuMP.@variable(model, 0.0 <= control_F[N, O, SF, G, T])
    JuMP.@variable(model, 0.0 <= release_location[N, T] <= 1.0)

    model.obj_dict[:control_limit_schedule_upper_F] =
        JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(undef, N, O, SF, G, T)
    model.obj_dict[:control_limit_schedule_lower_F] =
        JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(undef, N, O, SF, G, T)
    model.obj_dict[:control_limit_total_F] =
        JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(undef, N, O, SF, G, T)

    model.obj_dict[:control_limit_schedule_upper_M] =
        JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(undef, N, O, SM, G, T)
    model.obj_dict[:control_limit_schedule_lower_M] =
        JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(undef, N, O, SM, G, T)
    model.obj_dict[:control_limit_total_M] =
        JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(undef, N, O, SM, G, T)

    model.obj_dict[:control_equivalence] =
        JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(undef, N, O, T)

    # WARMSTART VARIABLES_1: Lifestages
    ###########################################
    stages = [Egg, Larva, Pupa, Male, Female]
    initialcond_dict = Dict(
        s => Dict(
            n => Dict{Int, Any}(
                c => Dict{Int, Matrix}(o => Matrix{Float64}(undef, 1, 1) for o in O) for c in C
            ) for n in N
        ) for s in stages
    )
    for (ix, node_name) in enumerate(N)
        for (cx, scenario) in enumerate(C)
            for (jx, organism) in enumerate(O)
                initialcond_dict[Egg][node_name][scenario][organism] =
                    initial_condition.x[ix].x[jx][SE, :]
                initialcond_dict[Larva][node_name][scenario][organism] =
                    initial_condition.x[ix].x[jx][(nE + 1):(nL + nE), :]
                initialcond_dict[Pupa][node_name][scenario][organism] =
                    initial_condition.x[ix].x[jx][(nE + nL + 1):(nE + nL + nP), :]
                initialcond_dict[Male][node_name][scenario][organism] =
                    initial_condition.x[ix].x[jx][nE + nL + nP + 1, :]'
                initialcond_dict[Female][node_name][scenario][organism] =
                    initial_condition.x[ix].x[jx][(nE + nL + nP + 2):end, :]
            end
        end
    end
    for node_name in N, scenario in C, organism in O, t in T #G[o]?
       # for genotype in G#[organism]

        JuMP.set_start_value.(
            E[node_name, scenario, organism, :, :, t].data, # no longer broadcasting over matrix, rather dict -> have to write as a loop (for x, set_start_value()) because now sparse 
            # make genes a matrix not a dict 
            initialcond_dict[Egg][node_name][scenario][organism],
        )
        JuMP.set_start_value.(
            L[node_name, scenario, organism, :, :, t].data,
            initialcond_dict[Larva][node_name][scenario][organism],
        )
        JuMP.set_start_value.(
            P[node_name, scenario, organism, :, :, t].data,
            initialcond_dict[Pupa][node_name][scenario][organism],
        )
        JuMP.set_start_value.(
            M[node_name, scenario, organism, :, :, t].data,
            initialcond_dict[Male][node_name][scenario][organism],
        )
        JuMP.set_start_value.(
            F[node_name, scenario, organism, :, :, t].data,
            initialcond_dict[Female][node_name][scenario][organism],
        )
       # end 
    end

    # EXPRESSIONS: Migration
    ###########################################
   # @show [G for o in O]
    #for species in node_species 
    # = get_migration(network, node_species)
    JuMP.@expression(
        model,
        migration_E[n in N, c in C, o in O, s in SE, g in G, t in T], # add  or brackets [G]?
        get_migration(network, node_species[o])[SE_map[s], g][n, :]' * E[:, c, o, s, g, t]
    )
    JuMP.@expression(
        model,
        migration_L[n in N, c in C, o in O, s in SL, g in G, t in T],
        get_migration(network, node_species[o])[SL_map[s], g][n, :]' * L[:, c, o, s, g, t]
    )
    JuMP.@expression(
        model,
        migration_P[n in N, c in C, o in O, s in SP, g in G, t in T],
        get_migration(network, node_species[o])[SP_map[s], g][n, :]' * P[:, c, o, s, g, t]
    )
    JuMP.@expression(
        model,
        migration_M[n in N, c in C, o in O, s in SM, g in G, t in T],
        get_migration(network, node_species[o])[SM_map[s], g][n, :]' * M[:, c, o, s, g, t]
    )
    # CONSTRAINTS_A: Life Stages
    ###########################################

    #### EGGS
    JuMP.@constraint(
        model,
        E_con_A0[n in N, c in C, o in O, s in [SE[1]], g in G, t in [T[1]]], # add brackets [G]?
        E[n, c, o, s, g, t] ==
        initialcond_dict[Egg][n][c][o][SE_map[s], g] + sum(
            (data[n]["scenario"][c]["organism"][o]["genetics"].likelihood[
                :,
                :,
                g,
            ] .* data[n]["scenario"][c]["organism"][o]["genetics"].Τ[
                :,
                :,
                g,
            ] .* data[n]["scenario"][c]["organism"][o]["genetics"].S[g] * data[n]["scenario"][c]["organism"][o]["genetics"].Β[g] .* F[
                n,
                c,
                o,
                :,
                :,
                t,
            ].data)[
                :,
                :,
            ],
        ) -
        E[n, c, o, s, g, t] * (
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Egg][n][c][o][g][t][1] *
            compute_density(densE, sum(E[n, c, o, :, :, t])) +
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Egg][n][c][o][g][t][2] *
            nE
        ) + migration_E[n, c, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        E_con_A1[n in N, c in C, o in O, s in [SE[1]], g in G, t in T[2:end]],
        E[n, c, o, s, g, t] ==
        E[n, c, o, s, g, t - 1] + sum(
            (data[n]["scenario"][c]["organism"][o]["genetics"].likelihood[
                :,
                :,
                g,
            ] .* data[n]["scenario"][c]["organism"][o]["genetics"].Τ[
                :,
                :,
                g,
            ] .* data[n]["scenario"][c]["organism"][o]["genetics"].S[g] * data[n]["scenario"][c]["organism"][o]["genetics"].Β[g] .* F[
                n,
                c,
                o,
                :,
                :,
                t,
            ].data)[
                :,
                :,
            ],
        ) -
        E[n, c, o, s, g, t] * (
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Egg][n][c][o][g][t][1] *
            compute_density(densE, sum(E[n, c, o, :, :, t])) +
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Egg][n][c][o][g][t][2] *
            nE
        ) + migration_E[n, c, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        E_con_B0[n in N, c in C, o in O, s in SE[2:end], g in G, t in [T[1]]],
        E[n, c, o, s, g, t] ==
        initial_condition.x[n].x[o][SE_map[s], g] +
        data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Egg][n][c][o][g][t][2] *
        nE *
        E[n, c, o, s - 1, g, t] -
        E[n, c, o, s, g, t] * (
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Egg][n][c][o][g][t][1] *
            compute_density(densE, sum(E[n, c, o, :, :, t])) +
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Egg][n][c][o][g][t][2] *
            nE
        ) + migration_E[n, c, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        E_con_B1[n in N, c in C, o in O, s in SE[2:end], g in G, t in T[2:end]],
        E[n, c, o, s, g, t] ==
        E[n, c, o, s, g, t - 1] +
        data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Egg][n][c][o][g][t][2] *
        nE *
        E[n, c, o, s - 1, g, t] -
        E[n, c, o, s, g, t] * (
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Egg][n][c][o][g][t][1] *
            compute_density(densE, sum(E[n, c, o, :, :, t])) +
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Egg][n][c][o][g][t][2] *
            nE
        ) + migration_E[n, c, o, s, g, t]
    )

    #### LARVAE
    JuMP.@constraint(
        model,
        L_con_A0[n in N, c in C, o in O, s in [SL[1]], g in G, t in [T[1]]],
        L[n, c, o, s, g, t] ==
        initial_condition.x[n].x[o][SL_map[s], g] +
        data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Egg][n][c][o][g][t][2] *
        nE *
        E[n, c, o, end, g, t] -
        L[n, c, o, s, g, t] * (
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Larva][n][c][o][g][t][1] *
            compute_density(densL, sum(L[n, c, o, :, :, t])) +
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Larva][n][c][o][g][t][2] *
            nL
        ) + migration_L[n, c, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        L_con_A1[n in N, c in C, o in O, s in [SL[1]], g in G, t in T[2:end]],
        L[n, c, o, s, g, t] ==
        L[n, c, o, s, g, t - 1] +
        data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Egg][n][c][o][g][t][2] *
        nE *
        E[n, c, o, end, g, t] -
        L[n, c, o, s, g, t] * (
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Larva][n][c][o][g][t][1] *
            compute_density(densL, sum(L[n, c, o, :, :, t])) +
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Larva][n][c][o][g][t][2] *
            nL
        ) + migration_L[n, c, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        L_con_B0[n in N, c in C, o in O, s in SL[2:end], g in G, t in [T[1]]],
        L[n, c, o, s, g, t] ==
        initial_condition.x[n].x[o][SL_map[s], g] +
        data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Larva][n][c][o][g][t][2] *
        nL *
        L[n, c, o, s - 1, g, t] -
        L[n, c, o, s, g, t] * (
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Larva][n][c][o][g][t][1] *
            compute_density(densL, sum(L[n, c, o, :, :, t])) +
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Larva][n][c][o][g][t][2] *
            nL
        ) + migration_L[n, c, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        L_con_B1[n in N, c in C, o in O, s in SL[2:end], g in G, t in T[2:end]],
        L[n, c, o, s, g, t] ==
        L[n, c, o, s, g, t - 1] +
        data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Larva][n][c][o][g][t][2] *
        nL *
        L[n, c, o, s - 1, g, t] -
        L[n, c, o, s, g, t] * (
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Larva][n][c][o][g][t][1] *
            compute_density(densL, sum(L[n, c, o, :, :, t])) +
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Larva][n][c][o][g][t][2] *
            nL
        ) + migration_L[n, c, o, s, g, t]
    )

    #### PUPAE
    JuMP.@constraint(
        model,
        P_con_A0[n in N, c in C, o in O, s in [SP[1]], g in G, t in [T[1]]],
        P[n, c, o, s, g, t] ==
        initial_condition.x[n].x[o][SP_map[s], g] +
        data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Larva][n][c][o][g][t][2] *
        nL *
        L[n, c, o, end, g, t] -
        P[n, c, o, s, g, t] * (
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][1] *
            compute_density(densP, sum(P[n, c, o, :, :, t])) +
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][2] *
            nP
        ) + migration_P[n, c, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        P_con_A1[n in N, c in C, o in O, s in [SP[1]], g in G, t in T[2:end]],
        P[n, c, o, s, g, t] ==
        P[n, c, o, s, g, t - 1] +
        data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Larva][n][c][o][g][t][2] *
        nL *
        L[n, c, o, end, g, t] -
        P[n, c, o, s, g, t] * (
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][1] *
            compute_density(densP, sum(P[n, c, o, :, :, t])) +
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][2] *
            nP
        ) + migration_P[n, c, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        P_con_B0[n in N, c in C, o in O, s in SP[2:end], g in G, t in [T[1]]],
        P[n, c, o, s, g, t] ==
        initial_condition.x[n].x[o][SP_map[s], g] +
        data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][2] *
        nP *
        P[n, c, o, s - 1, g, t] -
        P[n, c, o, s, g, t] * (
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][1] *
            compute_density(densP, sum(P[n, c, o, :, :, t])) +
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][2] *
            nP
        ) + migration_P[n, c, o, s, g, t]
    )

    JuMP.@constraint(
        model,
        P_con_B1[n in N, c in C, o in O, s in SP[2:end], g in G, t in T[2:end]],
        P[n, c, o, s, g, t] ==
        P[n, c, o, s, g, t - 1] +
        data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][2] *
        nP *
        P[n, c, o, s - 1, g, t] -
        P[n, c, o, s, g, t] * (
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][1] *
            compute_density(densP, sum(P[n, c, o, :, :, t])) +
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][2] *
            nP
        ) + migration_P[n, c, o, s, g, t]
    )

    #### MALES
    if slack_small || slack_large
        JuMP.@constraint(
            model,
            M_con_0[n in N, c in C, o in O, s in SM, g in G, t in [T[1]]],
            M[n, c, o, s, g, t] ==
            initial_condition.x[n].x[o][SM_map[s], g] +
            (1 - data[n]["scenario"][c]["organism"][o]["genetics"].Φ[g]) *
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][2] *
            nP *
            (P[n, c, o, end, g, t] + P_slack_positive[n, c, o, g, t]) * # - P_slack_negative[n, c, o, g, t])*
            data[n]["scenario"][c]["organism"][o]["genetics"].Ξ_m[g] -
            (1 + data[n]["scenario"][c]["organism"][o]["genetics"].Ω[g]) *
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Male][n][c][o][g][t][1] *
            M[n, c, o, s, g, t] *
            compute_density(densM, sum(M[n, c, o, :, :, t])) +
            migration_M[n, c, o, s, g, t]
        )
    else
        JuMP.@constraint(
            model,
            M_con_0[n in N, c in C, o in O, s in SM, g in G, t in [T[1]]],
            M[n, c, o, s, g, t] ==
            initial_condition.x[n].x[o][SM_map[s], g] +
            (1 - data[n]["scenario"][c]["organism"][o]["genetics"].Φ[g]) *
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][2] *
            nP *
            P[n, c, o, end, g, t] *
            data[n]["scenario"][c]["organism"][o]["genetics"].Ξ_m[g] -
            (1 + data[n]["scenario"][c]["organism"][o]["genetics"].Ω[g]) *
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Male][n][c][o][g][t][1] *
            M[n, c, o, s, g, t] *
            compute_density(densM, sum(M[n, c, o, :, :, t])) +
            migration_M[n, c, o, s, g, t]
        )
    end

    if slack_large
        JuMP.@constraint(
            model,
            M_con_1[n in N, c in C, o in O, s in SM, g in G, t in T[2:end]],
            M[n, c, o, s, g, t] ==
            M[n, c, o, s, g, t - 1] +
            (1 - data[n]["scenario"][c]["organism"][o]["genetics"].Φ[g]) *
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][2] *
            nP *
            (P[n, c, o, end, g, t] + P_slack_positive[n, c, o, g, t]) * #- P_slack_negative[n, c, o, g, t]) *
            data[n]["scenario"][c]["organism"][o]["genetics"].Ξ_m[g] -
            (1 + data[n]["scenario"][c]["organism"][o]["genetics"].Ω[g]) *
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Male][n][c][o][g][t][1] *
            M[n, c, o, s, g, t] *
            compute_density(densM, sum(M[n, c, o, :, :, t])) +
            control_M[n, o, s, g, t] +
            migration_M[n, c, o, s, g, t]
        )
    else
        JuMP.@constraint(
            model,
            M_con_1[n in N, c in C, o in O, s in SM, g in G, t in T[2:end]],
            M[n, c, o, s, g, t] ==
            M[n, c, o, s, g, t - 1] +
            (1 - data[n]["scenario"][c]["organism"][o]["genetics"].Φ[g]) *
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][2] *
            nP *
            P[n, c, o, end, g, t] *
            data[n]["scenario"][c]["organism"][o]["genetics"].Ξ_m[g] -
            (1 + data[n]["scenario"][c]["organism"][o]["genetics"].Ω[g]) *
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Male][n][c][o][g][t][1] *
            M[n, c, o, s, g, t] *
            compute_density(densM, sum(M[n, c, o, :, :, t])) +
            control_M[n, o, s, g, t] +
            migration_M[n, c, o, s, g, t]
        )
    end

    #### MATING
    JuMP.@constraint(
        model,
        mate_bound[n in N, c in C, o in O, g in G, t in T],
        M[n, c, o, 1, g, t] * data[n]["scenario"][c]["organism"][o]["genetics"].Η[g] <=
        (sum(
            M[n, c, o, 1, i, t] * data[n]["scenario"][c]["organism"][o]["genetics"].Η[i] for
            i in G
        ))
    )

    #### FEMALES:
    JuMP.@NLconstraint(
        model,
        F_con_0[n in N, c in C, o in O, s in SF, g in G, t in [T[1]]],
        F[n, c, o, s, g, t] ==
        initial_condition.x[n].x[o][SF_map[s], g] +
        (
            M[n, c, o, 1, g, t] * data[n]["scenario"][c]["organism"][o]["genetics"].Η[g] /
            (
                1e-6 + sum(
                    M[n, c, o, 1, i, t] *
                    data[n]["scenario"][c]["organism"][o]["genetics"].Η[i] for i in G
                )
            )
        ) * (
            data[n]["scenario"][c]["organism"][o]["genetics"].Φ[s] *
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][2] *
            nP *
            P[n, c, o, end, s, t] *
            data[n]["scenario"][c]["organism"][o]["genetics"].Ξ_f[s]
        ) -
        (1 + data[n]["scenario"][c]["organism"][o]["genetics"].Ω[g]) *
        data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Female][n][c][o][g][t][1] *
        F[n, c, o, s, g, t] + sum(get_migration(network, node_species[o])[SF_map[s], g][n, i] * F[i, c, o, s, g, t] for i in N)
    )

    JuMP.@NLconstraint(
        model,
        F_con_1[n in N, c in C, o in O, s in SF, g in G, t in T[2:end]],
        F[n, c, o, s, g, t] ==
        F[n, c, o, s, g, t - 1] +
        (
            M[n, c, o, 1, g, t] * data[n]["scenario"][c]["organism"][o]["genetics"].Η[g] /
            (
                1e-6 + sum(
                    M[n, c, o, 1, i, t] *
                    data[n]["scenario"][c]["organism"][o]["genetics"].Η[i] for i in G
                )
            )
        ) * (
            data[n]["scenario"][c]["organism"][o]["genetics"].Φ[s] *
            data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Pupa][n][c][o][g][t][2] *
            nP *
            P[n, c, o, end, s, t] *
            data[n]["scenario"][c]["organism"][o]["genetics"].Ξ_f[s]
        ) -
        (1 + data[n]["scenario"][c]["organism"][o]["genetics"].Ω[g]) *
        data[n]["scenario"][c]["organism"][o]["stage_temperature_response"][Female][n][c][o][g][t][1] *
        F[n, c, o, s, g, t] +
        control_F[n, o, s, g, t] +
        sum(get_migration(network, node_species[o])[SF_map[s], g][n, i] * F[i, c, o, s, g, t] for i in N)
    )

    # CONSTRAINTS_B: Controls
    ###########################################
    for org_key in keys(organism_data["organism"])
        homozygous_modified = organism_data["organism"][org_key]["homozygous_modified"]
        wildtype = organism_data["organism"][org_key]["wildtype"]

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

    end 
    
    return model
end

"""
    create_decision_model(node::Node, tspan; node_strategy::DataStructures.OrderedDict, species::Type{<:Species}=AedesAegypti,do_binary::Bool=false, optimizer=nothing)

Build mathematical program. Problem created as an NLP (do_binary=false) or MINLP (do_binary=true). NB: `Node` is recreated as a `Network` object internally; this does not change the problem but is relevant for data exploration as it adds one index layer to the formatted results.
"""
function create_decision_model(node::Node, tspan; kwargs...)
    node_name = get_name(node)
    network = Network(node_name, node)
    return create_decision_model(network, tspan; kwargs...)
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
    # Re-set constraints by fixing controls to zero
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
