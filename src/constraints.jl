
function _optimization_info(network::Network, tspan::Tuple)
    optinfo_dict = Dict()

    # Eco/bio data
    ####################
    spatialtempresponse_dict = Dict()

    optinfo_dict["time_horizon"] = tspan
    optinfo_dict["network_name"] = get_name(network)
    optinfo_dict["total_node_count"] = count_nodes(network)
    optinfo_dict["network_total_organism_count"] = 0

    for (ix, node) in enumerate(values(get_nodes(network)))
        optinfo_dict[ix] = Dict{String, Any}()
        optinfo_dict[ix]["organism"] = Dict()
        optinfo_dict[ix]["temperature"] = get_temperature(node).values

        for (jx, organism) in enumerate(get_organisms(node))
            species_count = count_organisms(node)
            optinfo_dict["network_total_organism_count"] += species_count
            optinfo_dict[ix]["node_organism_count"] = species_count
            optinfo_dict[ix]["organism"][jx] = Dict{String, Any}()
            optinfo_dict[ix]["organism"][jx]["genetics"] = get_genetics(node, organism)
            optinfo_dict[ix]["organism"][jx]["gene_count"] = count_genotypes(node, organism)
            optinfo_dict[ix]["organism"][jx]["homozygous_modified"] =
                get_homozygous_modified(node, organism)
            optinfo_dict[ix]["organism"][jx]["wildtype"] = get_wildtype(node, organism)
            optinfo_dict[ix]["organism"][jx]["total_stage_count"] = 0
            optinfo_dict[ix]["organism"][jx]["substage_count"] = Dict()
            optinfo_dict[ix]["organism"][jx]["stage_density"] = Dict()

            for (key_lifestage, life_stage) in get_lifestages(node, organism)
                substage_count = count_substages(node, organism, key_lifestage)
                optinfo_dict[ix]["organism"][jx]["substage_count"][key_lifestage] =
                    substage_count
                optinfo_dict[ix]["organism"][jx]["total_stage_count"] += substage_count
                optinfo_dict[ix]["organism"][jx]["stage_temperature_response"] = Dict()
                optinfo_dict[ix]["organism"][jx]["stage_density"][key_lifestage] =
                    get_density(node, organism, key_lifestage)

                for lx in 1:count_genotypes(network, get_name(node), organism)
                    stage_dict = get!(spatialtempresponse_dict, key_lifestage, Dict())
                    node_dict = get!(stage_dict, ix, Dict())
                    org_response_dict = get!(node_dict, jx, Dict())
                    org_response_dict[lx] = [
                        temperature_effect(temp, life_stage) for
                        temp in node.temperature.values
                    ]
                    optinfo_dict[ix]["organism"][jx]["stage_temperature_response"] =
                        spatialtempresponse_dict
                end
            end
        end
    end
    return optinfo_dict
end

"""
    mutable struct ReleaseStrategy
        release_this_gene_index::Union{Nothing, Int64}=nothing
        release_this_life_stage::Union{Nothing, Type{<:LifeStage}, Type{Female}}=nothing
        release_location_force::Union{Nothing, Bool}=nothing
        release_start_time::Union{Nothing,Int64}=nothing
        release_end_time::Union{Nothing,Int64}=nothing
        release_time_interval::Int64=1
        release_size_max_per_timestep::Union{Int64,Float64}=9e9
        release_max_over_timehorizon::Union{Int64,Float64}=9e9

    end

Data defining the operational constraints for each node.

# Arguments

  - `release_this_gene_index`: Gene index to be released. Generally defined by `get_homozygous_modified`.
  - `release_this_life_stage`: Lifestage to be released. Varies according to genetic technology.
  - `release_location_force`: Specify locations where releases are obligatory; only applicable when decision model is being run as an MINLP.
  - `release_start_time`: Timestep (day) that releases are permitted to start.
  - `release_end_time`: Timestep (day) that releases are required to end.
  - `release_time_interval`: Minimum timestep interval (in days) permitted for releases.
  - `release_size_max_per_timestep`: Maximum number of organisms that may be released on a single day.
  - `release_max_over_timehorizon`: Maximum number of organisms that may be released over the problem horizon.
"""
Base.@kwdef mutable struct ReleaseStrategy
    release_this_gene_index::Union{Nothing, Int64} = nothing
    release_this_life_stage::Union{Nothing, Type{<:LifeStage}, Type{Female}} = nothing
    release_location_force::Union{Nothing, Bool} = nothing
    release_start_time::Union{Nothing, Int64} = nothing
    release_end_time::Union{Nothing, Int64} = nothing
    release_time_interval::Int64 = 1
    release_size_max_per_timestep::Union{Int64, Float64} = 9e9
    release_max_over_timehorizon::Union{Int64, Float64} = 9e9
end

function _add_constraint(
    model::JuMP.Model,
    optinfo_dict::Dict,
    life_stage::Type{<:Male},
    node_strategy::Dict,
    homozygous_modified,
    wildtype,
)
    data = optinfo_dict["constraints"]
    max_overall = node_strategy[1].release_max_over_timehorizon #TODO: fix
    control_M = model[:control_M]
    sets = model.obj_dict[:Sets]
    N = sets[:N]
    O = sets[:O]
    SM = sets[:SM]
    G = sets[:G]
    T = sets[:T]

    # SCHEDULE
    JuMP.@constraint(
        model,
        control_limit_schedule[n in N, o in O, s in SM, g in G, t in T],
        control_M[n, o, 1, g, t] <= data[n]["org_con"][o]["stage_con"][life_stage][g][t]
    )

    # MAXIMUM
    JuMP.@constraint(model, control_limit_total, sum(control_M) <= max_overall)

    # FIX ALTERNATIVE(S) 
    control_F = model[:control_F]
    JuMP.fix.(control_F, 0.0; force=true)
    JuMP.fix.(control_F[:, :, wildtype, :, :], 0.0; force=true)
    JuMP.fix.(control_F[:, :, :, homozygous_modified, :], 0.0; force=true)

    return
end

function _add_constraint(
    model::JuMP.Model,
    optinfo_dict::Dict,
    life_stage::Type{<:Female},
    node_strategy::Dict,
    homozygous_modified,
    wildtype,
)
    data = optinfo_dict["constraints"]
    max_overall = node_strategy[1].release_max_over_timehorizon #TODO: fix
    control_F = model[:control_F]
    sets = model.obj_dict[:Sets]
    N = sets[:N]
    O = sets[:O]
    SF = sets[:SF]
    G = sets[:G]
    T = sets[:T]

    # SCHEDULE
    JuMP.@constraint(
        model,
        control_limit_schedule[n in N, o in O, s in SF, g in G, t in T],
        control_F[n, o, wildtype, homozygous_modified, t] <=
        data[n]["org_con"][o]["stage_con"][life_stage][g][t]
    )

    # MAXIMUM
    JuMP.@constraint(model, control_limit_total, sum(control_F) <= max_overall)

    # FIX ALTERNATIVE(S)
    control_M = model[:control_M]
    JuMP.fix.(control_M, 0.0; force=true)

    return
end

function _calculate_release_constraints(
    network::Network,
    tspan::Tuple,
    homozygous_modified,
    wildtype,
    do_binary::Bool,
    model::JuMP.Model,
    optinfo_dict::Dict,
    node_strategy::Dict{Int64, ReleaseStrategy},
)
    releaseconstraints_dict = Dict()

    nodes = get_nodes(network)
    for (ix, node) in enumerate(values(nodes))
        releaseconstraints_dict[ix] = Dict{String, Any}()
        releaseconstraints_dict[ix]["org_con"] = Dict()

        for (jx, organism) in enumerate(get_organisms(node))
            releaseconstraints_dict[ix]["org_con"][jx] = Dict{String, Any}()
            releaseconstraints_dict[ix]["org_con"][jx]["stage_con"] = Dict()

            for (key_lifestage, life_stage) in get_lifestages(node, organism)
                releaseconstraints_dict[ix]["org_con"][jx]["stage_con"][key_lifestage] =
                    Dict()
                for lx in 1:count_genotypes(network, get_name(node), organism)
                    releaseconstraints_dict[ix]["org_con"][jx]["stage_con"][key_lifestage][lx] =
                        zeros(length(tspan[1]:tspan[2]))
                end
            end
        end
    end

    # Checks
    ####################
    for (node_number, strategy) in node_strategy
        if strategy.release_this_gene_index === nothing
            @info "Please define the index of the genotype to be released."
        end

        if strategy.release_this_life_stage === nothing
            @info "Please define the lifestage of the organism to be released."
        end

        if strategy.release_location_force !== nothing && do_binary == false
            @warn "Release locations cannot be specified because `do_binary` is set to $(do_binary)."
        end

        if strategy.release_end_time !== nothing && strategy.release_end_time > tspan[2]
            @warn "Release end time specified in node $(entry[1]) occurs outside the time horizon (tspan = $(tspan))."
        end

        if all(
            isequal(
                node_strategy[1].release_max_over_timehorizon,
                strategy.release_max_over_timehorizon,
            ),
        ) == false
            @warn "Entries for the `ReleaseStrategy` field `release_max_over_timehorizon` do not match; this field should be equivalent for all nodes in the network."
        end

        strategy.release_start_time === nothing ? (strategy.release_start_time = tspan[1]) :
        strategy.release_start_time
        strategy.release_end_time === nothing ? (strategy.release_end_time = tspan[2]) :
        strategy.release_end_time
    end

    for (node_number, strategy) in node_strategy

        # TODO: fix
        for tx in
            (strategy.release_start_time):(strategy.release_time_interval):(strategy.release_end_time)
            releaseconstraints_dict[node_number]["org_con"][1]["stage_con"][strategy.release_this_life_stage][strategy.release_this_gene_index][tx] =
                strategy.release_size_max_per_timestep
        end
    end

    optinfo_dict["constraints"] = releaseconstraints_dict

    _add_constraint(
        model,
        optinfo_dict,
        node_strategy[1].release_this_life_stage,
        node_strategy,
        homozygous_modified,
        wildtype,
    )

    return optinfo_dict
end
