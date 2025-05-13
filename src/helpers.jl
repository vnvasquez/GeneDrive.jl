########################################
#         Solve Dynamic Model          #
########################################
"""
    solve_dynamic_model(node::Node, algorithm, tspan)

Return ODE model solution for single node problem with no releases.
"""
function solve_dynamic_model(node::Node, algorithm, tspan)
    @info("No releases specified.")

    tstops = Vector()
    callbacks = Vector()
    collected_callback_set = []

    if isa(node.temperature, TimeSeriesTemperature)
        tempseries = [
            TemperatureSeriesData(
                node,
                collect(tspan[1]:tspan[2]),
                node.temperature.values,
            ),
        ]
        for temp in tempseries
            push!(tstops, temp.times...)
            push!(callbacks, temp.set.discrete_callbacks...)
        end
        collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))

    elseif isa(node.temperature, ScenarioTemperature) # && selected_scenario !== nothing
        tempseries = [
            TemperatureSeriesData(
                node,
                collect(tspan[1]:tspan[2]),
                node.temperature.values[:, node.temperature.selected_scenario],
            ),
        ]

        for temp in tempseries
            push!(tstops, temp.times...)
            push!(callbacks, temp.set.discrete_callbacks...)
        end
        collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))

    else
        collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))
    end

    inputs = ExogenousInputs(node)
    u0, dens = init_node!(node)
    @info("Simulation initialized. Beginning model run with $(typeof(node.temperature)).")

    prob = diffeq.ODEProblem(population_model_node, u0, tspan, (node, inputs))
    sol = diffeq.solve(
        prob,
        algorithm,
        callback=collected_callback_set,
        tstops=unique!(tstops),
        reltol=1e-9,
    )
    @info(" Model run complete.")

    return sol
end

"""
    solve_dynamic_model(node::Node, releases::Union{Vector{Release},Vector{ProportionalRelease}}, algorithm, tspan)

Return ODE model solution for single node problem with releases.
"""
function solve_dynamic_model(
    node::Node,
    releases::Union{Vector{Release}, Vector{ProportionalRelease}},
    algorithm,
    tspan,
)
    tstops = Vector()
    callbacks = Vector()
    collected_callback_set = []

    if isa(node.temperature, TimeSeriesTemperature)
        tempseries = [
            TemperatureSeriesData(
                node,
                collect(tspan[1]:tspan[2]),
                node.temperature.values,
            ),
        ]
        for release in releases
            push!(tstops, release.times...)
            push!(callbacks, release.callbacks...)
        end
        for temp in tempseries
            push!(tstops, temp.times...)
            push!(callbacks, temp.set.discrete_callbacks...)
        end
        collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))

    elseif isa(node.temperature, ScenarioTemperature) # && selected_scenario !== nothing
        tempseries = [
            TemperatureSeriesData(
                node,
                collect(tspan[1]:tspan[2]),
                node.temperature.values[:, node.temperature.selected_scenario],
            ),
        ]

        for release in releases
            push!(tstops, release.times...)
            push!(callbacks, release.callbacks...)
        end
        for temp in tempseries
            push!(tstops, temp.times...)
            push!(callbacks, temp.set.discrete_callbacks...)
        end
        collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))

    else
        for release in releases
            push!(tstops, release.times...)
            push!(callbacks, release.callbacks...)
        end
        collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))
    end
    inputs = ExogenousInputs(node)

    u0, dens = init_node!(node)
    @info(
        "Simulation initialized. Beginning model run with releases and $(typeof(node.temperature))."
    )
    prob = diffeq.ODEProblem(population_model_node, u0, tspan, (node, inputs))

    sol = diffeq.solve(
        prob,
        algorithm,
        callback=collected_callback_set,
        tstops=unique!(tstops),
        reltol=1e-9,
    )
    @info(" Model run complete.")

    return sol
end

"""
    solve_dynamic_model(node::Node, shocks::TemperatureShockData, algorithm, tspan)

Return ODE model solution for single node problem with temperature shocks.
"""
function solve_dynamic_model(node::Node, shocks::TemperatureShockData, algorithm, tspan)
    @info("No releases specified.")

    tstops = Vector{Float64}()
    callbacks = Vector()
    collected_callback_set = []

    if isa(node.temperature, TimeSeriesTemperature)
        shock_data = zeros(length(node.temperature.values))
        for (ix, t) in enumerate(tspan[1]:tspan[2])
            for (jx, t_range) in enumerate(shocks.times)
                if t_range[1] <= t <= t_range[2]
                    shock_data[ix] = shocks.values[jx]
                end
            end
        end
        tempseries = [
            TemperatureSeriesData(
                node,
                collect(tspan[1]:tspan[2]),
                node.temperature.values + shock_data,
            ),
        ]
        for temp in tempseries
            push!(tstops, temp.times...)
            push!(callbacks, temp.set.discrete_callbacks...)
        end
        collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))
    else
        for times in shocks.times
            push!(tstops, times...)
        end
        push!(callbacks, shocks.set.discrete_callbacks...)
    end
    collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))
    inputs = ExogenousInputs(node)
    u0, dens = init_node!(node)
    @info(
        "Simulation initialized. Beginning model run with shocks and $(typeof(node.temperature))."
    )

    prob = diffeq.ODEProblem(population_model_node, u0, tspan, (node, inputs))
    sol = diffeq.solve(
        prob,
        algorithm,
        callback=collected_callback_set,
        tstops=unique!(tstops),
        reltol=1e-9,
        dt=1.0,
    )
    @info(" Model run complete.")

    return sol
end

"""
    solve_dynamic_model(node::Node, releases::Union{Vector{Release},Vector{ProportionalRelease}}, shocks::TemperatureShockData, algorithm, tspan)

Return ODE model solution for single node problem with releases and temperature shocks.
"""
function solve_dynamic_model(
    node::Node,
    releases::Union{Vector{Release}, Vector{ProportionalRelease}},
    shocks::TemperatureShockData,
    algorithm,
    tspan,
)
    tstops = Vector{Float64}()
    callbacks = Vector()
    collected_callback_set = []

    if isa(node.temperature, TimeSeriesTemperature)
        shock_data = zeros(length(node.temperature.values))
        for (ix, t) in enumerate(tspan[1]:tspan[2])
            for (jx, t_range) in enumerate(shocks.times)
                if t_range[1] <= t <= t_range[2]
                    shock_data[ix] = shocks.values[jx]
                end
            end
        end
        tempseries = [
            TemperatureSeriesData(
                node,
                collect(tspan[1]:tspan[2]),
                node.temperature.values + shock_data,
            ),
        ]
        for temp in tempseries
            push!(tstops, temp.times...)
            push!(callbacks, temp.set.discrete_callbacks...)
        end
        for release in releases
            push!(tstops, release.times...)
            push!(callbacks, release.callbacks...)
        end
        collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))
    else
        for release in releases
            push!(tstops, release.times...)
            push!(callbacks, release.callbacks...)
        end
        for times in shocks.times
            push!(tstops, times...)
        end
        push!(callbacks, shocks.set.discrete_callbacks...)
        collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))
    end
    inputs = ExogenousInputs(node)
    u0, dens = init_node!(node)
    @info(
        "Simulation initialized. Beginning model run with releases, shocks, and $(typeof(node.temperature))."
    )

    prob = diffeq.ODEProblem(population_model_node, u0, tspan, (node, inputs))
    sol = diffeq.solve(
        prob,
        algorithm,
        callback=collected_callback_set,
        tstops=unique!(tstops),
        reltol=1e-9,
        dt=1.0,
    )
    @info(" Model run complete.")

    return sol
end

"""
    solve_dynamic_model(network::Network, algorithm, tspan)

Return ODE model solution for network problem with no releases.
"""
function solve_dynamic_model(network::Network, algorithm, tspan)
    @info("No releases specified for any node of the network $(get_name(network)).")

    tstops = Vector()
    callbacks = Vector()
    collected_callback_set = []

    for (key, node) in network.nodes
        if isa(node.temperature, TimeSeriesTemperature)
            tempseries = [
                TemperatureSeriesData(
                    node,
                    collect(tspan[1]:tspan[2]),
                    node.temperature.values,
                ),
            ]

            for temp in tempseries
                push!(tstops, temp.times...)
                push!(callbacks, temp.set.discrete_callbacks...)
            end
            collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))

        elseif isa(node.temperature, ScenarioTemperature)
            tempseries = [
                TemperatureSeriesData(
                    node,
                    collect(tspan[1]:tspan[2]),
                    node.temperature.values[:, node.temperature.selected_scenario],
                ),
            ]
            for temp in tempseries
                push!(tstops, temp.times...)
                push!(callbacks, temp.set.discrete_callbacks...)
            end
            collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))

        else
            collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))
        end
    end
    inputs = ExogenousInputs(network)
    u0_net, dens_net = init_network!(network)
    @info(
        "Simulation initialized. Beginning model run in $(count_nodes(network))-node network $(get_name(network))."
    )

    prob = diffeq.ODEProblem(population_model_network, u0_net, tspan, (network, inputs))
    sol_net = diffeq.solve(
        prob,
        algorithm,
        callback=collected_callback_set,
        tstops=unique!(tstops),
        reltol=1e-9,
    )
    @info(" Model run complete.")

    return sol_net
end

"""
    solve_dynamic_model(network::Network, releases::Union{Vector{Release}, Vector{ProportionalRelease}}, algorithm, tspan)

Return ODE model solution for network problem with releases.
"""
function solve_dynamic_model(
    network::Network,
    releases::Union{Vector{Release}, Vector{ProportionalRelease}},
    algorithm,
    tspan,
)
    @info("Releases specified for the following nodes in the network: ")
    for release in releases
        println("$(release.node)")
    end

    tstops = Vector()
    callbacks = Vector()
    collected_callback_set = []

    for (key, node) in network.nodes
        if isa(node.temperature, TimeSeriesTemperature)
            tempseries = [
                TemperatureSeriesData(
                    node,
                    collect(tspan[1]:tspan[2]),
                    node.temperature.values,
                ),
            ]

            for temp in tempseries
                push!(tstops, temp.times...)
                push!(callbacks, temp.set.discrete_callbacks...)
            end
            for release in releases
                push!(tstops, release.times...)
                push!(callbacks, release.callbacks...)
            end
            collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))

        elseif isa(node.temperature, ScenarioTemperature)
            tempseries = [
                TemperatureSeriesData(
                    node,
                    collect(tspan[1]:tspan[2]),
                    node.temperature.values[:, node.temperature.selected_scenario],
                ),
            ]
            for temp in tempseries
                push!(tstops, temp.times...)
                push!(callbacks, temp.set.discrete_callbacks...)
            end
            for release in releases
                push!(tstops, release.times...)
                push!(callbacks, release.callbacks...)
            end
            collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))

        else
            for release in releases
                push!(tstops, release.times...)
                push!(callbacks, release.callbacks...)
            end
            collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))
        end
    end
    inputs = ExogenousInputs(network)
    u0_net, dens_net = init_network!(network)
    @info(
        "Simulation initialized. Beginning model run (with releases) in $(count_nodes(network))-node network $(get_name(network))."
    )

    prob = diffeq.ODEProblem(population_model_network, u0_net, tspan, (network, inputs))
    sol_net = diffeq.solve(
        prob,
        algorithm,
        callback=collected_callback_set,
        tstops=unique!(tstops),
        reltol=1e-9,
    )
    @info(" Model run complete.")

    return sol_net
end

"""
    solve_dynamic_model(network::Network, shocks::Vector{TemperatureShockData}, algorithm, tspan)

Return ODE model solution for network problem with temperature shocks.
"""
function solve_dynamic_model(
    network::Network,
    shocks::Vector{TemperatureShockData},
    algorithm,
    tspan,
)
    @info("No releases specified for any node of the network $(get_name(network)).")

    @info("Temperature shocks specified for the following nodes in the network: ")
    for shock in shocks
        println("$(shock.node)")
    end

    tstops = Vector()
    callbacks = Vector()
    collected_callback_set = []

    for (key, node) in network.nodes
        if isa(node.temperature, TimeSeriesTemperature)
            shock_data = zeros(length(node.temperature.values))
            for (ix, t) in enumerate(tspan[1]:tspan[2])
                for shock in shocks
                    for (jx, t_range) in enumerate(shock.times)
                        if t_range[1] <= t <= t_range[2]
                            shock_data[ix] = shock.values[jx]
                        end
                    end
                end
            end
            tempseries = [
                TemperatureSeriesData(
                    node,
                    collect(tspan[1]:tspan[2]),
                    node.temperature.values,
                ),
            ]
            for temp in tempseries
                push!(tstops, temp.times...)
                push!(callbacks, temp.set.discrete_callbacks...)
            end
            collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))
        else
            for shock in shocks
                for times in shock.times
                    push!(tstops, times...)
                end
                push!(callbacks, shock.set.discrete_callbacks...)
            end
            collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))
        end
    end
    inputs = ExogenousInputs(network)
    u0_net, dens_net = init_network!(network)
    @info(
        "Simulation initialized. Beginning model run (with temperature shocks) in $(count_nodes(network))-node network $(get_name(network))."
    )

    prob = diffeq.ODEProblem(population_model_network, u0_net, tspan, (network, inputs))
    sol_net = diffeq.solve(
        prob,
        algorithm,
        callback=collected_callback_set,
        tstops=unique!(tstops),
        reltol=1e-9,
    )
    @info(" Model run complete.")

    return sol_net
end

"""
    solve_dynamic_model(network::Network, releases::Union{Vector{Release},Vector{ProportionalRelease}}, shocks::Vector{TemperatureShockData}, algorithm, tspan)

Return ODE model solution for network problem with releases and temperature shocks.
"""
function solve_dynamic_model(
    network::Network,
    releases::Union{Vector{Release}, Vector{ProportionalRelease}},
    shocks::Vector{TemperatureShockData},
    algorithm,
    tspan,
)
    @info("Releases specified for the following nodes in the network: ")
    for release in releases
        println("$(release.node)")
    end

    @info("Temperature shocks specified for the following nodes in the network: ")
    for shock in shocks
        println("$(shock.node)")
    end

    tstops = Vector()
    callbacks = Vector()
    collected_callback_set = []

    for (key, node) in network.nodes
        if isa(node.temperature, TimeSeriesTemperature)
            shock_data = zeros(length(node.temperature.values))
            for (ix, t) in enumerate(tspan[1]:tspan[2])
                for shock in shocks
                    for (jx, t_range) in enumerate(shock.times)
                        if t_range[1] <= t <= t_range[2]
                            shock_data[ix] = shock.values[jx]
                        end
                    end
                end
            end
            tempseries = [
                TemperatureSeriesData(
                    node,
                    collect(tspan[1]:tspan[2]),
                    node.temperature.values,
                ),
            ]
            for temp in tempseries
                push!(tstops, temp.times...)
                push!(callbacks, temp.set.discrete_callbacks...)
            end
            for release in releases
                push!(tstops, release.times...)
                push!(callbacks, release.callbacks...)
            end
            collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))
        else
            for shock in shocks
                for times in shock.times
                    push!(tstops, times...)
                end
                push!(callbacks, shock.set.discrete_callbacks...)
            end
            for release in releases
                push!(tstops, release.times...)
                push!(callbacks, release.callbacks...)
            end
            collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))
        end
    end
    inputs = ExogenousInputs(network)
    u0_net, dens_net = init_network!(network)
    @info(
        "Simulation initialized. Beginning model run (with releases and temperature shocks) in $(count_nodes(network))-node network $(get_name(network))."
    )

    prob = diffeq.ODEProblem(population_model_network, u0_net, tspan, (network, inputs))
    sol_net = diffeq.solve(
        prob,
        algorithm,
        callback=collected_callback_set,
        tstops=unique!(tstops),
        reltol=1e-9,
    )
    @info(" Model run complete.")

    return sol_net
end

"""
    solve_scenarios_dynamic_model(node::Node, algorithm, tspan, scenarios_of_interest::Vector{Int})

Return ODE model solution for node problem with scenarios and no releases.
"""
function solve_scenarios_dynamic_model(
    node::Node,
    algorithm,
    tspan,
    scenarios_of_interest::Vector{Int},
)
    solutions = []
    for scenario in scenarios_of_interest
        set_scenario!(node, scenario)
        @info "Solving temperature scenario #$(node.temperature.selected_scenario)..."
        sol = solve_dynamic_model(node, algorithm, tspan)
        push!(solutions, sol)
    end

    return solutions
end

"""
    solve_scenarios_dynamic_model(node::Node, releases, algorithm, tspan, scenarios_of_interest::Vector{Int})

Return ODE model solution for node problem with scenarios and releases.
"""
function solve_scenarios_dynamic_model(
    node::Node,
    releases::Union{Vector{Release}, Vector{ProportionalRelease}},
    algorithm,
    tspan,
    scenarios_of_interest::Vector{Int},
)
    solutions = []
    for scenario in scenarios_of_interest
        set_scenario!(node, scenario)
        @info "Solving temperature scenario #$(node.temperature.selected_scenario)..."
        sol = solve_dynamic_model(node, releases, algorithm, tspan)
        push!(solutions, sol)
    end

    return solutions
end

"""
    solve_scenarios_dynamic_model(network::Network, algorithm, tspan, scenarios_of_interest::Vector{Int})

Return ODE model solution for network problem with scenarions and no releases.
"""
function solve_scenarios_dynamic_model(
    network::Network,
    algorithm,
    tspan,
    scenarios_of_interest::Vector{Int},
)
    solutions = []
    for scenario in scenarios_of_interest
        for (key, node) in network.nodes
            set_scenario!(node, scenario)
        end
        @info "Solving temperature scenario #$(scenario)..."
        sol = solve_dynamic_model(network, algorithm, tspan)
        push!(solutions, sol)
    end

    return solutions
end

"""
    solve_scenarios_dynamic_model(network::Network, releases, algorithm, tspan, scenarios_of_interest::Vector{Int})

Return ODE model solution for node problem with scenarios and releases.
"""
function solve_scenarios_dynamic_model(
    network::Network,
    releases::Union{Vector{Release}, Vector{ProportionalRelease}},
    algorithm,
    tspan,
    scenarios_of_interest::Vector{Int},
)
    solutions = []
    for scenario in scenarios_of_interest
        for (key, node) in network.nodes
            set_scenario!(node, scenario)
        end
        @info "Solving temperature scenario #$(scenario)..."
        sol = solve_dynamic_model(network, releases, algorithm, tspan)
        push!(solutions, sol)
    end

    return solutions
end

########################################
#        Format Dynamic Results        #
########################################

function _count_substages(node::Node, species::Type{<:Species}, stage::Type{<:LifeStage})
    return node.organisms[species].all_stages[stage].n
end

function _count_substages(node::Node, species::Type{<:Species}, stage::Type{Female})
    substages = node.organisms[species].all_stages[stage].n
    gN = count_genotypes(node, species)

    return substages * gN
end

"""
    format_dynamic_model_results(node::Node, sol)

Return dictionary containing ODE model results. Indexed per organism, life stage, and genotype.
"""
function format_dynamic_model_results(node::Node, sol)
    node_results = sol.prob.p[1]

    organisms = get_organisms(node)
    gN = count_genotypes(node_results, get_organisms(node_results)[1]) #TODO: FIX 
    timesteps = length(sol.t)

    results = Dict{String, Dict{String, Dict{String, Matrix{Float64}}}}()

    for (index, org) in enumerate(organisms)
        org_key = "$(org)"
        results[org_key] = Dict{String, Dict{String, Matrix{Float64}}}()

        stages = get_lifestages(node, org)
        genotypes = get_genotypes(node, org)

        for (index, stage) in enumerate(keys(stages))
            stage_key = "$(stage)"
            sN = _count_substages(node, org, stage)
            results[org_key][stage_key] = Dict{String, Matrix{Float64}}()

            for gene in genotypes
                gene_key = "$(gene.genotype)"
                results[org_key][stage_key][gene_key] = zeros(Float64, sN, timesteps)
            end
        end
    end

    timesteps = length(sol.t)
    for t in 1:(timesteps - 1)
        for (index, org) in enumerate(organisms)
            sN_ = 1
            org_index = index
            org_key = "$(org)"
            mat = sol.u[t].x[org_index]

            genotypes = get_genotypes(node, org)
            stages = get_lifestages(node, org)

            for (index, stage) in enumerate(keys(stages))
                stage_index = index
                stage_key = "$(stage)"
                sN = _count_substages(node, org, stage)

                for (index, gene) in enumerate(genotypes)
                    gene_index = index
                    gene_key = "$(gene.genotype)"
                    range = sN_:(sN_ + sN - 1)

                    if stage_key != "Female"
                        values = mat[range, gene_index]
                    else
                        values = transpose(mat[range, gene_index])
                    end

                    results[org_key][stage_key][gene_key][1:sN, t] = values
                end

                sN_ += sN
            end
        end
    end
    return results
end

"""
    format_dynamic_model_results(network::Network, sol)

Return dictionary containing ODE model results. Indexed per node, organism, life stage, and genotype.
"""
function format_dynamic_model_results(network::Network, sol)
    timesteps = length(sol.t)
    results = Dict{String, Dict{String, Dict{String, Dict{String, Matrix{Float64}}}}}()

    for (index, node) in enumerate(values(network.nodes))
        node_key = "$(get_name(node))"
        results[node_key] = Dict{String, Dict{String, Dict{String, Matrix{Float64}}}}()
        organisms = get_organisms(node)
        gN = count_genotypes(node, get_organisms(node)[1])

        for (index, org) in enumerate(organisms)
            org_key = "$(org)"
            results[node_key][org_key] = Dict{String, Dict{String, Matrix{Float64}}}()

            stages = get_lifestages(node, org)
            genotypes = get_genotypes(node, org)

            for (index, stage) in enumerate(keys(stages))
                stage_key = "$(stage)"
                sN = _count_substages(node, org, stage)
                results[node_key][org_key][stage_key] = Dict{String, Matrix{Float64}}()

                for gene in genotypes
                    gene_key = "$(gene.genotype)"
                    results[node_key][org_key][stage_key][gene_key] =
                        zeros(Float64, sN, timesteps)
                end
            end
        end
    end

    for t in 1:(timesteps - 1)
        for (index, node) in enumerate(values(network.nodes))
            node_index = index
            node_key = "$(get_name(node))"
            organisms = get_organisms(node)

            for (index, org) in enumerate(organisms)#enumerate(keys(organisms))
                sN_ = 1
                org_index = index
                org_key = "$(org)"
                mat = sol.u[t].x[node_index].x[org_index]

                genotypes = get_genotypes(node, org)
                stages = get_lifestages(node, org)

                for (index, stage) in enumerate(keys(stages))
                    stage_index = index
                    stage_key = "$(stage)"
                    sN = _count_substages(node, org, stage)

                    for (index, gene) in enumerate(genotypes)
                        gene_index = index
                        gene_key = "$(gene.genotype)"
                        range = sN_:(sN_ + sN - 1)

                        if stage_key != "Female"
                            values = mat[range, gene_index]
                        else
                            values = transpose(mat[range, gene_index])
                        end

                        results[node_key][org_key][stage_key][gene_key][1:sN, t] = values
                    end

                    sN_ += sN
                end
            end
        end
    end
    return results
end

########################################
#        Format Decision Results       #
########################################
"""
    format_decision_model_results(sol)

Return dictionary containing optimization model results. Indexed per node, organism, life stage, and genotype.
"""
function format_decision_model_results(sol)
    results_dict = Dict()
    results_dict[:Time] = collect(sol.obj_dict[:Sets][:T])
    for (var_key, var_val) in sol.obj_dict
        if any(occursin("release_location", String(var_key)))
            @info("Excluding $(var_key) variable from results.")
            continue
        end

        if eltype(var_val) == JuMP.VariableRef
            for n in axes(var_val)[1]
                if length(axes(var_val)) == 6 # bc scenarios
                    for c in axes(var_val)[2]
                        for o in axes(var_val)[3]
                            df = DataFrames.DataFrame()
                            for g in axes(var_val[n, c, o, :, :, :])[2]
                                col_symbol = Symbol(var_key, "_G$(g)")
                                if occursin("F", String(var_key))
                                    df[!, col_symbol] = sum(
                                        JuMP.value.(var_val[n, c, o, g, :, :]).data,
                                        dims=1,
                                    )[
                                        1,
                                        :,
                                    ]
                                    key_symbol = Symbol(
                                        "node_$(n)_scenario_$(c)_organism_$(o)_$(var_key)",
                                    )
                                    results_dict[key_symbol] = df
                                    continue
                                end
                                df[!, col_symbol] =
                                    sum(JuMP.value.(var_val[n, c, o, :, g, :]).data, dims=1)[
                                        1,
                                        :,
                                    ]
                                key_symbol = Symbol(
                                    "node_$(n)_scenario_$(c)_organism_$(o)_$(var_key)",
                                )
                                results_dict[key_symbol] = df
                            end
                        end
                    end

                elseif length(axes(var_val)) == 5 && occursin("slack", String(var_key))
                    for c in axes(var_val)[2]
                        for o in axes(var_val)[3]
                            df = DataFrames.DataFrame()
                            for g in axes(var_val[n, c, o, :, :])[2]
                                col_symbol = Symbol(var_key, "_G$(g)")
                                df[!, col_symbol] =
                                    sum(JuMP.value.(var_val[n, c, o, g, :]).data, dims=1)[
                                        1,
                                        :,
                                    ]
                                key_symbol = Symbol(
                                    "node_$(n)_scenario_$(c)_organism_$(o)_$(var_key)",
                                )
                                results_dict[key_symbol] = df
                            end
                        end
                    end

                    #elseif length(axes(var_val)) == 5 && occursin("control", String(var_key))
                else#if occursin("control", String(var_key))
                    for o in axes(var_val)[2]
                        df = DataFrames.DataFrame()
                        for g in axes(var_val[n, o, :, :, :])[2]
                            col_symbol = Symbol(var_key, "_G$(g)")
                            if occursin("F", String(var_key))
                                df[!, col_symbol] =
                                    sum(JuMP.value.(var_val[n, o, g, :, :]).data, dims=1)[
                                        1,
                                        :,
                                    ]
                                key_symbol = Symbol("node_$(n)_organism_$(o)_$(var_key)")
                                results_dict[key_symbol] = df
                                continue
                            end
                            df[!, col_symbol] =
                                sum(JuMP.value.(var_val[n, o, :, g, :]).data, dims=1)[1, :]
                            key_symbol = Symbol("node_$(n)_organism_$(o)_$(var_key)")
                            results_dict[key_symbol] = df
                        end
                    end
                end
            end

            if any(occursin("release_location", String(var_key)))
                # return release_locations .= release_location
            end
        end
    end
    return results_dict
end

########################################
#     Plot Select Dynamic Results      #
########################################

"""
    plot_dynamic_mendelian_females(node::Node, sol)

Return visualization of adult female population dynamics across all genotypes.
"""
function plot_dynamic_mendelian_females(node::Node, sol)
    results = format_dynamic_model_results(node, sol)

    mendelian_base_F_1 = []
    mendelian_base_F_2 = []
    mendelian_base_F_3 = []

    for k in keys(node.organisms)
        k = string(k)

        mendelian_base_F_1 =
            results[k]["Female"]["AA"][1, :] .+ results[k]["Female"]["Aa"][1, :] .+
            results[k]["Female"]["aa"][1, :]

        mendelian_base_F_2 =
            results[k]["Female"]["AA"][2, :] .+ results[k]["Female"]["Aa"][2, :] .+
            results[k]["Female"]["aa"][2, :]

        mendelian_base_F_3 =
            results[k]["Female"]["AA"][3, :] .+ results[k]["Female"]["Aa"][3, :] .+
            results[k]["Female"]["aa"][3, :]
    end

    timesteps = sol.t[1:(end - 1)]

    traces = PlotlyJS.GenericTrace{Dict{Symbol, Any}}[]
    push!(
        traces,
        PlotlyJS.scatter(;
            name="AA",
            x=timesteps,
            y=mendelian_base_F_1,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
        ),
    )
    push!(
        traces,
        PlotlyJS.scatter(;
            name="Aa",
            x=timesteps,
            y=mendelian_base_F_2[1:(end - 1)],
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dashdot",
        ),
    )

    push!(
        traces,
        PlotlyJS.scatter(;
            name="aa",
            x=timesteps,
            y=mendelian_base_F_3[1:(end - 1)],
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dot",
        ),
    )
    p = PlotlyJS.plot(
        traces,
        PlotlyJS.Layout(
            title=PlotlyJS.attr(text="Females", x=0.5, xanchor="center"),
            xaxis_title="Time [Days]",
            yaxis_title="Population [Count]",
            width=800,
            height=450,
            font_size=12,
            fillcolor="transparent",
            xaxis=PlotlyJS.attr(tickfont_size=14),
            yaxis=PlotlyJS.attr(tickfont_size=14),
            legend=PlotlyJS.attr(
                yanchor="bottom",
                y=-0.3,
                xanchor="center",
                x=0.5,
                orientation="h",
                font_size=11.2,
            ),
        ),
    )
end

"""
    plot_dynamic_wolbachia_females(node::Node, sol)

Return visualization of adult female population dynamics across all genotypes.
"""
function plot_dynamic_wolbachia_females(node::Node, sol)
    results = format_dynamic_model_results(node, sol)

    wolbachia_base_F_1 = []
    wolbachia_base_F_2 = []

    for k in keys(node.organisms)
        k = string(k)

        wolbachia_base_F_1 =
            results[k]["Female"]["WW"][1, :] .+ results[k]["Female"]["ww"][1, :]
        wolbachia_base_F_2 =
            results[k]["Female"]["WW"][2, :] .+ results[k]["Female"]["ww"][2, :]
    end

    timesteps = sol.t[1:(end - 1)]

    traces = PlotlyJS.GenericTrace{Dict{Symbol, Any}}[]
    push!(
        traces,
        PlotlyJS.scatter(;
            name="Wildtype",
            x=timesteps,
            y=wolbachia_base_F_2,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
        ),
    )
    push!(
        traces,
        PlotlyJS.scatter(;
            name="Infected",
            x=timesteps,
            y=wolbachia_base_F_1[1:(end - 1)],
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dashdot",
        ),
    )
    p = PlotlyJS.plot(
        traces,
        PlotlyJS.Layout(
            title=PlotlyJS.attr(text="Females", x=0.5, xanchor="center"),
            xaxis_title="Time [Days]",
            yaxis_title="Population [Count]",
            width=800,
            height=450,
            font_size=12,
            fillcolor="transparent",
            xaxis=PlotlyJS.attr(tickfont_size=14),
            yaxis=PlotlyJS.attr(tickfont_size=14),
            legend=PlotlyJS.attr(
                yanchor="bottom",
                y=-0.3,
                xanchor="center",
                x=0.5,
                orientation="h",
                font_size=11.2,
            ),
        ),
    )
end

"""
    plot_dynamic_ridl_females(node::Node, sol)

Return visualization of adult female population dynamics across all genotypes.
"""
function plot_dynamic_ridl_females(node::Node, sol)
    results = format_dynamic_model_results(node, sol)

    ridl_base_F_1 = []
    ridl_base_F_2 = []
    ridl_base_F_3 = []

    for k in keys(node.organisms)
        k = string(k)

        ridl_base_F_1 =
            results[k]["Female"]["WW"][1, :] .+ results[k]["Female"]["WR"][1, :] .+
            results[k]["Female"]["RR"][1, :]

        ridl_base_F_2 =
            results[k]["Female"]["WW"][2, :] .+ results[k]["Female"]["WR"][2, :] .+
            results[k]["Female"]["RR"][2, :]

        ridl_base_F_3 =
            results[k]["Female"]["WW"][3, :] .+ results[k]["Female"]["WR"][3, :] .+
            results[k]["Female"]["RR"][3, :]
    end

    timesteps = sol.t[1:(end - 1)]

    traces = PlotlyJS.GenericTrace{Dict{Symbol, Any}}[]
    push!(
        traces,
        PlotlyJS.scatter(;
            name="WW",
            x=timesteps,
            y=ridl_base_F_1,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
        ),
    )
    push!(
        traces,
        PlotlyJS.scatter(;
            name="WR",
            x=timesteps,
            y=ridl_base_F_2[1:(end - 1)],
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dashdot",
        ),
    )

    push!(
        traces,
        PlotlyJS.scatter(;
            name="RR",
            x=timesteps,
            y=ridl_base_F_3[1:(end - 1)],
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dot",
        ),
    )
    p = PlotlyJS.plot(
        traces,
        PlotlyJS.Layout(
            title=PlotlyJS.attr(text="Females", x=0.5, xanchor="center"),
            xaxis_title="Time [Days]",
            yaxis_title="Population [Count]",
            width=800,
            height=450,
            font_size=12,
            fillcolor="transparent",
            xaxis=PlotlyJS.attr(tickfont_size=14),
            yaxis=PlotlyJS.attr(tickfont_size=14),
            legend=PlotlyJS.attr(
                yanchor="bottom",
                y=-0.3,
                xanchor="center",
                x=0.5,
                orientation="h",
                font_size=11.2,
            ),
        ),
    )
end

"""
    plot_dynamic_mcr_females(node::Node, sol)

Return visualization of adult female population dynamics across all genotypes.
"""
function plot_dynamic_mcr_females(node::Node, sol)
    results = format_dynamic_model_results(node, sol)

    mcr_base_F_1 = []
    mcr_base_F_2 = []
    mcr_base_F_3 = []
    mcr_base_F_4 = []
    mcr_base_F_5 = []
    mcr_base_F_6 = []

    for k in keys(node.organisms)
        k = string(k)

        mcr_base_F_1 =
            results[k]["Female"]["HH"][1, :] .+ results[k]["Female"]["Hh"][1, :] .+
            results[k]["Female"]["HR"][1, :] .+ results[k]["Female"]["hh"][1, :] .+
            results[k]["Female"]["hR"][1, :] .+ results[k]["Female"]["RR"][1, :]

        mcr_base_F_2 =
            results[k]["Female"]["HH"][2, :] .+ results[k]["Female"]["Hh"][2, :] .+
            results[k]["Female"]["HR"][2, :] .+ results[k]["Female"]["hh"][2, :] .+
            results[k]["Female"]["hR"][2, :] .+ results[k]["Female"]["RR"][2, :]

        mcr_base_F_3 =
            results[k]["Female"]["HH"][3, :] .+ results[k]["Female"]["Hh"][3, :] .+
            results[k]["Female"]["HR"][3, :] .+ results[k]["Female"]["hh"][3, :] .+
            results[k]["Female"]["hR"][3, :] .+ results[k]["Female"]["RR"][3, :]

        mcr_base_F_4 =
            results[k]["Female"]["HH"][4, :] .+ results[k]["Female"]["Hh"][4, :] .+
            results[k]["Female"]["HR"][4, :] .+ results[k]["Female"]["hh"][4, :] .+
            results[k]["Female"]["hR"][4, :] .+ results[k]["Female"]["RR"][4, :]

        mcr_base_F_5 =
            results[k]["Female"]["HH"][5, :] .+ results[k]["Female"]["Hh"][5, :] .+
            results[k]["Female"]["HR"][5, :] .+ results[k]["Female"]["hh"][5, :] .+
            results[k]["Female"]["hR"][5, :] .+ results[k]["Female"]["RR"][5, :]

        mcr_base_F_6 =
            results[k]["Female"]["HH"][6, :] .+ results[k]["Female"]["Hh"][6, :] .+
            results[k]["Female"]["HR"][6, :] .+ results[k]["Female"]["hh"][6, :] .+
            results[k]["Female"]["hR"][6, :] .+ results[k]["Female"]["RR"][6, :]
    end

    timesteps = sol.t[1:(end - 1)]

    traces = PlotlyJS.GenericTrace{Dict{Symbol, Any}}[]
    push!(
        traces,
        PlotlyJS.scatter(;
            name="HH",
            x=timesteps,
            y=mcr_base_F_1,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
        ),
    )
    push!(
        traces,
        PlotlyJS.scatter(;
            name="Hh",
            x=timesteps,
            y=mcr_base_F_2[1:(end - 1)],
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dashdot",
        ),
    )

    push!(
        traces,
        PlotlyJS.scatter(;
            name="HR",
            x=timesteps,
            y=mcr_base_F_3[1:(end - 1)],
            mode="lines",
            line_shape="linear",
            line_color="red",
            fillcolor="transparent",
            line_width=2,
            line_dash="dashdot",
        ),
    )

    push!(
        traces,
        PlotlyJS.scatter(;
            name="hh",
            x=timesteps,
            y=mcr_base_F_4[1:(end - 1)],
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dot",
        ),
    )

    push!(
        traces,
        PlotlyJS.scatter(;
            name="hR",
            x=timesteps,
            y=mcr_base_F_5[1:(end - 1)],
            mode="lines",
            line_shape="linear",
            line_color="red",
            fillcolor="transparent",
            line_width=2,
            line_dash="dot",
        ),
    )

    push!(
        traces,
        PlotlyJS.scatter(;
            name="RR",
            x=timesteps,
            y=mcr_base_F_6[1:(end - 1)],
            mode="lines",
            line_shape="linear",
            line_color="red",
            fillcolor="transparent",
            line_width=2,
        ),
    )

    p = PlotlyJS.plot(
        traces,
        PlotlyJS.Layout(
            title=PlotlyJS.attr(text="Females", x=0.5, xanchor="center"),
            xaxis_title="Time [Days]",
            yaxis_title="Population [Count]",
            width=800,
            height=450,
            font_size=12,
            fillcolor="transparent",
            xaxis=PlotlyJS.attr(tickfont_size=14),
            yaxis=PlotlyJS.attr(tickfont_size=14),
            legend=PlotlyJS.attr(
                yanchor="bottom",
                y=-0.3,
                xanchor="center",
                x=0.5,
                orientation="h",
                font_size=11.2,
            ),
        ),
    )
end

########################################
#     Plot Select Decision Results      #
########################################

"""
    plot_decision_mendelian_females(sol)

Return visualization of adult female population dynamics across all genotypes.
"""
function plot_decision_mendelian_females(sol)
    results = format_decision_model_results(sol)

    df = results[:node_1_organism_1_F]
    timesteps = results[:Time]

    traces = PlotlyJS.GenericTrace{Dict{Symbol, Any}}[]
    push!(
        traces,
        PlotlyJS.scatter(;
            name="AA",
            x=timesteps,
            y=df.F_G1,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
        ),
    )
    push!(
        traces,
        PlotlyJS.scatter(;
            name="Aa",
            x=timesteps,
            y=df.F_G2,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dashdot",
        ),
    )

    push!(
        traces,
        PlotlyJS.scatter(;
            name="aa",
            x=timesteps,
            y=df.F_G3,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dot",
        ),
    )
    p = PlotlyJS.plot(
        traces,
        PlotlyJS.Layout(
            title=PlotlyJS.attr(text="Females", x=0.5, xanchor="center"),
            xaxis_title="Time [Days]",
            yaxis_title="Population [Count]",
            width=800,
            height=450,
            font_size=12,
            fillcolor="transparent",
            xaxis=PlotlyJS.attr(tickfont_size=14),
            yaxis=PlotlyJS.attr(tickfont_size=14),
            legend=PlotlyJS.attr(
                yanchor="bottom",
                y=-0.3,
                xanchor="center",
                x=0.5,
                orientation="h",
                font_size=11.2,
            ),
        ),
    )
end

"""
    plot_decision_ridl_females(sol)

Return visualization of adult female population dynamics across all genotypes.
"""
function plot_decision_ridl_females(sol)
    results = format_decision_model_results(sol)

    df = results[:node_1_organism_1_F]
    timesteps = results[:Time]

    traces = PlotlyJS.GenericTrace{Dict{Symbol, Any}}[]
    push!(
        traces,
        PlotlyJS.scatter(;
            name="WW",
            x=timesteps,
            y=df.F_G1,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
        ),
    )
    push!(
        traces,
        PlotlyJS.scatter(;
            name="WR",
            x=timesteps,
            y=df.F_G2,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dashdot",
        ),
    )

    push!(
        traces,
        PlotlyJS.scatter(;
            name="RR",
            x=timesteps,
            y=df.F_G3,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dot",
        ),
    )
    p = PlotlyJS.plot(
        traces,
        PlotlyJS.Layout(
            title=PlotlyJS.attr(text="Females", x=0.5, xanchor="center"),
            xaxis_title="Time [Days]",
            yaxis_title="Population [Count]",
            width=800,
            height=450,
            font_size=12,
            fillcolor="transparent",
            xaxis=PlotlyJS.attr(tickfont_size=14),
            yaxis=PlotlyJS.attr(tickfont_size=14),
            legend=PlotlyJS.attr(
                yanchor="bottom",
                y=-0.3,
                xanchor="center",
                x=0.5,
                orientation="h",
                font_size=11.2,
            ),
        ),
    )
end

"""
    plot_decision_mcr_females(sol)

Return visualization of adult female population dynamics across all genotypes.
"""
function plot_decision_mcr_females(sol)
    results = format_decision_model_results(sol)

    df = results[:node_1_organism_1_F]
    timesteps = results[:Time]

    traces = PlotlyJS.GenericTrace{Dict{Symbol, Any}}[]
    push!(
        traces,
        PlotlyJS.scatter(;
            name="HH",
            x=timesteps,
            y=df.F_G1,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
        ),
    )
    push!(
        traces,
        PlotlyJS.scatter(;
            name="Hh",
            x=timesteps,
            y=df.F_G2,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dashdot",
        ),
    )

    push!(
        traces,
        PlotlyJS.scatter(;
            name="HR",
            x=timesteps,
            y=df.F_G3,
            mode="lines",
            line_shape="linear",
            line_color="red",
            fillcolor="transparent",
            line_width=2,
            line_dash="dashdot",
        ),
    )

    push!(
        traces,
        PlotlyJS.scatter(;
            name="hh",
            x=timesteps,
            y=df.F_G4,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dot",
        ),
    )

    push!(
        traces,
        PlotlyJS.scatter(;
            name="hR",
            x=timesteps,
            y=df.F_G5,
            mode="lines",
            line_shape="linear",
            line_color="red",
            fillcolor="transparent",
            line_width=2,
            line_dash="dot",
        ),
    )

    push!(
        traces,
        PlotlyJS.scatter(;
            name="RR",
            x=timesteps,
            y=df.F_G6,
            mode="lines",
            line_shape="linear",
            line_color="red",
            fillcolor="transparent",
            line_width=2,
        ),
    )
    p = PlotlyJS.plot(
        traces,
        PlotlyJS.Layout(
            title=PlotlyJS.attr(text="Females", x=0.5, xanchor="center"),
            xaxis_title="Time [Days]",
            yaxis_title="Population [Count]",
            width=800,
            height=450,
            font_size=12,
            fillcolor="transparent",
            xaxis=PlotlyJS.attr(tickfont_size=14),
            yaxis=PlotlyJS.attr(tickfont_size=14),
            legend=PlotlyJS.attr(
                yanchor="bottom",
                y=-0.3,
                xanchor="center",
                x=0.5,
                orientation="h",
                font_size=11.2,
            ),
        ),
    )
end

"""
    plot_decision_wolbachia_females(sol)

Return visualization of adult female population dynamics across all genotypes.
"""
function plot_decision_wolbachia_females(sol)
    results = format_decision_model_results(sol)

    df = results[:node_1_organism_1_F]
    timesteps = results[:Time]

    traces = PlotlyJS.GenericTrace{Dict{Symbol, Any}}[]
    push!(
        traces,
        PlotlyJS.scatter(;
            name="WW",
            x=timesteps,
            y=df.F_G1,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
        ),
    )
    push!(
        traces,
        PlotlyJS.scatter(;
            name="ww",
            x=timesteps,
            y=df.F_G2,
            mode="lines",
            line_shape="linear",
            line_color="green",
            fillcolor="transparent",
            line_width=2,
            line_dash="dashdot",
        ),
    )

    p = PlotlyJS.plot(
        traces,
        PlotlyJS.Layout(
            title=PlotlyJS.attr(text="Females", x=0.5, xanchor="center"),
            xaxis_title="Time [Days]",
            yaxis_title="Population [Count]",
            width=800,
            height=450,
            font_size=12,
            fillcolor="transparent",
            xaxis=PlotlyJS.attr(tickfont_size=14),
            yaxis=PlotlyJS.attr(tickfont_size=14),
            legend=PlotlyJS.attr(
                yanchor="bottom",
                y=-0.3,
                xanchor="center",
                x=0.5,
                orientation="h",
                font_size=11.2,
            ),
        ),
    )
end
