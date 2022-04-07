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

    if isa(node.temperature, TimeSeriesTemperature) == true
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
    solve_dynamic_model(node::Node, releases::Vector, algorithm, tspan)

Return ODE model solution for single node problem with releases.
"""
function solve_dynamic_model(node::Node, releases::Vector, algorithm, tspan)
    tstops = Vector()
    callbacks = Vector()
    collected_callback_set = []

    if isa(node.temperature, TimeSeriesTemperature) == true
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

    if isa(node.temperature, TimeSeriesTemperature) == true
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
    solve_dynamic_model(node::Node, releases::Vector, shocks::TemperatureShockData, algorithm, tspan)

Return ODE model solution for single node problem with releases and temperature shocks.
"""
function solve_dynamic_model(
    node::Node,
    releases::Vector,
    shocks::TemperatureShockData,
    algorithm,
    tspan,
)
    tstops = Vector{Float64}()
    callbacks = Vector()
    collected_callback_set = []

    if isa(node.temperature, TimeSeriesTemperature) == true
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
        if isa(node.temperature, TimeSeriesTemperature) == true
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
    solve_dynamic_model(network::Network, releases::Vector, algorithm, tspan)

Return ODE model solution for network problem with releases.
"""
function solve_dynamic_model(network::Network, releases::Vector, algorithm, tspan)
    @info("Releases specified for the following nodes in the network: ")
    for release in releases
        println("$(release.node)")
    end

    tstops = Vector()
    callbacks = Vector()
    collected_callback_set = []

    for (key, node) in network.nodes
        if isa(node.temperature, TimeSeriesTemperature) == true
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
    solve_dynamic_model(network::Network, shocks::Vector, algorithm, tspan)

Return ODE model solution for network problem with temperature shocks.
"""
function solve_dynamic_model(network::Network, shocks::Vector, algorithm, tspan)
    @info("No releases specified for any node of the network $(get_name(network)).")

    @info("Temperature shocks specified for the following nodes in the network: ")
    for shock in shocks
        println("$(shock.node)")
    end

    tstops = Vector()
    callbacks = Vector()
    collected_callback_set = []

    for (key, node) in network.nodes
        if isa(node.temperature, TimeSeriesTemperature) == true
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
    solve_dynamic_model(network::Network, releases::Vector, shocks::Vector, algorithm, tspan)

Return ODE model solution for network problem with releases and temperature shocks.
"""
function solve_dynamic_model(
    network::Network,
    releases::Vector,
    shocks::Vector,
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
        if isa(node.temperature, TimeSeriesTemperature) == true
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
    gN = count_genotypes(node_results, get_organisms(node_results)[1])
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
    for (var_key, var_val) in sol.obj_dict
        if any(occursin("release_location", String(var_key)))
            @info("Excluding $(var_key) variable from results.")
            @info("FYI, $(var_key) is a $(var_val)")
            continue
        end

        if eltype(var_val) == VariableRef
            for n in axes(var_val)[1]
                for o in axes(var_val)[2]
                    df = DataFrame()
                    if length(axes(var_val)) == 5
                        for g in axes(var_val[n, o, :, :, :])[2]
                            col_symbol = Symbol(var_key, "_G$(g)")
                            if occursin("F", String(var_key))
                                df[!, col_symbol] =
                                    sum(value.(var_val[n, o, g, :, :]).data, dims=1)[1, :]
                                continue
                            end
                            df[!, col_symbol] =
                                sum(value.(var_val[n, o, :, g, :]).data, dims=1)[1, :]
                        end
                    end
                    if length(axes(var_val)) == 4
                        for g in axes(var_val)[3]
                            col_symbol = Symbol(var_key, "_G$(g)")
                        end
                    end
                    results_dict[Symbol("node_$(n)_organism_$(o)_$(var_key)")] = df
                    if any(occursin("release_location", String(var_key)))
                        return release_locations .= release_location
                    end
                end
            end
        end
    end
    return results_dict
end
