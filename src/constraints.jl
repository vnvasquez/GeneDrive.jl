
function _optimization_info(network::Network, tspan::Tuple)
    optinfo_dict = Dict()

    # Eco/bio data
    ####################
    spatialtempresponse_dict = Dict()

    optinfo_dict["time_horizon"] = tspan
    optinfo_dict["network_name"] = get_name(network)
    optinfo_dict["total_node_count"] = count_nodes(network)
    optinfo_dict["total_scenario_count"] = 0
    optinfo_dict["network_total_organism_count"] = 0

    for (ix, node) in enumerate(values(get_nodes(network)))
        optinfo_dict[ix] = Dict{String, Any}()
        optinfo_dict[ix]["scenario"] = Dict()

        for (cx, prob) in enumerate(get_probability(node.temperature))
            optinfo_dict[ix]["scenario"][cx] = Dict()
            optinfo_dict[ix]["scenario"][cx]["organism"] = Dict()

            temperature_vals = get_temperature_scenarios(node.temperature)[:, cx]
            optinfo_dict["total_scenario_count"] =
                count_temperature_scenarios(node.temperature)
            optinfo_dict[ix]["scenario"][cx]["temperature"] = temperature_vals
            optinfo_dict[ix]["scenario"][cx]["probability"] = prob

            for (jx, organism) in enumerate(get_organisms(node))
                species_count = count_organisms(node)
                optinfo_dict["network_total_organism_count"] += species_count
                optinfo_dict[ix]["scenario"][cx]["node_organism_count"] = species_count

                optinfo_dict[ix]["scenario"][cx]["organism"][jx] = Dict()
                optinfo_dict[ix]["scenario"][cx]["organism"][jx]["genetics"] =
                    get_genetics(node, organism)
                optinfo_dict[ix]["scenario"][cx]["organism"][jx]["gene_count"] =
                    count_genotypes(node, organism)
                optinfo_dict[ix]["scenario"][cx]["organism"][jx]["homozygous_modified"] =
                    get_homozygous_modified(node, organism)
                optinfo_dict[ix]["scenario"][cx]["organism"][jx]["wildtype"] =
                    get_wildtype(node, organism)
                optinfo_dict[ix]["scenario"][cx]["organism"][jx]["total_stage_count"] = 0
                optinfo_dict[ix]["scenario"][cx]["organism"][jx]["substage_count"] = Dict()
                optinfo_dict[ix]["scenario"][cx]["organism"][jx]["stage_density"] = Dict()

                for (key_lifestage, life_stage) in get_lifestages(node, organism)
                    substage_count = count_substages(node, organism, key_lifestage)
                    optinfo_dict[ix]["scenario"][cx]["organism"][jx]["substage_count"][key_lifestage] =
                        substage_count
                    optinfo_dict[ix]["scenario"][cx]["organism"][jx]["total_stage_count"] +=
                        substage_count
                    optinfo_dict[ix]["scenario"][cx]["organism"][jx]["stage_temperature_response"] =
                        Dict()
                    optinfo_dict[ix]["scenario"][cx]["organism"][jx]["stage_density"][key_lifestage] =
                        get_density(node, organism, key_lifestage)

                    #@show optinfo_dict[ix]["scenario"][cx]["organism"][jx]["stage_density"][key_lifestage]
                    for lx in 1:count_genotypes(network, get_name(node), organism)

                        stage_dict = get!(spatialtempresponse_dict, key_lifestage, Dict())
                        node_dict = get!(stage_dict, ix, Dict())
                        scenario_dict = get!(node_dict, cx, Dict())
                        org_response_dict = get!(scenario_dict, jx, Dict())

                        org_response_dict[lx] = [
                            temperature_effect(temp, life_stage) for
                            temp in temperature_vals
                        ]

                        optinfo_dict[ix]["scenario"][cx]["organism"][jx]["stage_temperature_response"] =
                            spatialtempresponse_dict
                    end
                end
            end
        end
    end
    return optinfo_dict
end

"""
    mutable struct ReleaseStrategy
        release_this_gene_index::Union{Nothing, Int64}=nothing
        release_this_life_stage=nothing
        release_location_force::Union{Nothing, Bool}=nothing
        release_start_time::Union{Nothing,Int64}=nothing
        release_end_time::Union{Nothing,Int64}=nothing
        release_time_interval::Int64=1
        release_size_min_per_timestep::Union{Int64, Float64}=0.0
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
  - `release_size_min_per_timestep`: Minimum number of organisms that may be released on a single day.
  - `release_size_max_per_timestep`: Maximum number of organisms that may be released on a single day.
  - `release_max_over_timehorizon`: Maximum number of organisms that may be released over the problem horizon.
"""
Base.@kwdef mutable struct ReleaseStrategy
    release_this_gene_index::Union{Nothing, Int64} = nothing
    release_this_life_stage = nothing
    release_location_force::Union{Nothing, Bool} = nothing
    release_start_time::Union{Nothing, Int64} = nothing
    release_end_time::Union{Nothing, Int64} = nothing
    release_time_interval::Int64 = 1
    release_size_min_per_timestep::Union{Int64, Float64} = 0.0
    release_size_max_per_timestep::Union{Int64, Float64} = 9e9
    release_max_over_timehorizon::Union{Int64, Float64} = 9e9
end

function _add_constraint(
    model::JuMP.Model,
    do_binary::Bool,
    optinfo_dict::Dict,
    node_number::Int64,
    life_stage::Type{<:Male},
    node_strategy::ReleaseStrategy,
    homozygous_modified,
    wildtype,
)
    data = optinfo_dict["constraints"]
    max_per_timestep = node_strategy.release_size_max_per_timestep
    min_per_timestep = node_strategy.release_size_min_per_timestep
    max_overall = node_strategy.release_max_over_timehorizon #TODO: fix
    control_M = model[:control_M]
    release_location = model[:release_location]
    sets = model.obj_dict[:Sets]
    N = sets[:N]
    O = sets[:O]
    SM = sets[:SM]
    G = sets[:G]
    T = sets[:T]

    CONTAINER_control_limit_schedule_upper_M =
        model.obj_dict[:control_limit_schedule_upper_M]
    CONTAINER_control_limit_schedule_lower_M =
        model.obj_dict[:control_limit_schedule_lower_M]
    CONTAINER_control_limit_total_M = model.obj_dict[:control_limit_total_M]

    range =
        (node_strategy.release_start_time):(node_strategy.release_time_interval):(node_strategy.release_end_time)

    # SCHEDULE
    for n in node_number, t in T
        if t ∈ range
            if do_binary
                #@info "made binary for node $n and time $t"
                JuMP.set_binary(release_location[n, t])
            end
            for o in O, s in SM, g in G
                max_v = g == homozygous_modified ? max_per_timestep : 0.0
                min_v = g == homozygous_modified ? min_per_timestep : 0.0
                #@info "max limit used in $s, $g is $max_v"
                CONTAINER_control_limit_schedule_upper_M[n, o, s, g, t] = JuMP.@constraint(
                    model,
                    # max
                    control_M[n, o, s, g, t] <= release_location[n, t] * max_v
                )
                CONTAINER_control_limit_schedule_lower_M[n, o, s, g, t] = JuMP.@constraint(
                    model,
                    # min
                    control_M[n, o, s, g, t] >= release_location[n, t] * min_v
                )
            end
        else
            #@info "fixed to 0.0 node $n and time $t"
            JuMP.fix(release_location[n, t], 0.0; force=true)
            for o in O, s in SM, g in G
                CONTAINER_control_limit_schedule_upper_M[n, o, s, g, t] =
                    JuMP.@constraint(model, control_M[n, o, s, g, t] <= 0)
                CONTAINER_control_limit_schedule_lower_M[n, o, s, g, t] =
                    JuMP.@constraint(model, control_M[n, o, s, g, t] >= 0)
            end
        end
    end

    # MAXIMUM
    CONTAINER_control_limit_total_M = JuMP.@constraint(model, sum(control_M) <= max_overall)

    control_F = model[:control_F]
    JuMP.fix.(control_F, 0.0; force=true)

    return
end

function _add_constraint(
    model::JuMP.Model,
    do_binary::Bool,
    optinfo_dict::Dict,
    node_number::Int64,
    life_stage::Type{<:Female},
    node_strategy::ReleaseStrategy,
    homozygous_modified,
    wildtype,
)
    data = optinfo_dict["constraints"]
    max_per_timestep = node_strategy.release_size_max_per_timestep
    min_per_timestep = node_strategy.release_size_min_per_timestep
    max_overall = node_strategy.release_max_over_timehorizon #TODO: fix
    control_F = model[:control_F]
    release_location = model[:release_location]
    sets = model.obj_dict[:Sets]
    N = sets[:N]
    O = sets[:O]
    SF = sets[:SF]
    G = sets[:G]
    T = sets[:T]

    CONTAINER_control_limit_schedule_upper_F =
        model.obj_dict[:control_limit_schedule_upper_F]
    CONTAINER_control_limit_schedule_lower_F =
        model.obj_dict[:control_limit_schedule_lower_F]
    CONTAINER_control_limit_total_F = model.obj_dict[:control_limit_total_F]

    range =
        (node_strategy.release_start_time):(node_strategy.release_time_interval):(node_strategy.release_end_time)

    # SCHEDULE
    for n in node_number, t in T
        if t ∈ range
            if do_binary
                #@info "made FEMALES binary for node $n and time $t"
                JuMP.set_binary(release_location[n, t])
            end
            for o in O, s in SF, g in G
                max_v = g == homozygous_modified ? max_per_timestep : 0.0
                min_v = g == homozygous_modified ? min_per_timestep : 0.0
                #@info "max limit used in $o, $s, $g is FEMALES = $max_v"
                CONTAINER_control_limit_schedule_upper_F[n, o, s, g, t] = JuMP.@constraint(
                    model,
                    # max
                    control_F[n, o, g, s, t] <= release_location[n, t] * max_v
                )
                CONTAINER_control_limit_schedule_lower_F[n, o, s, g, t] = JuMP.@constraint(
                    model,
                    # min
                    control_F[n, o, g, s, t] >= release_location[n, t] * min_v
                )
            end
        else
            #@info "fixed FEMALES to 0.0 node $n and time $t"
            JuMP.fix(release_location[n, t], 0.0; force=true)
            for o in O, s in SF, g in G
                CONTAINER_control_limit_schedule_upper_F[n, o, s, g, t] =
                    JuMP.@constraint(model, control_F[n, o, g, s, t] <= 0)
                CONTAINER_control_limit_schedule_lower_F[n, o, s, g, t] =
                    JuMP.@constraint(model, control_F[n, o, g, s, t] >= 0)
            end
        end
    end

    # MAXIMUM
    CONTAINER_control_limit_total_F = JuMP.@constraint(
        model,
        sum(control_F) <= max_overall * MALE_FEMALE_RELEASE_FRACTION
    )

    control_M = model[:control_M]
    JuMP.fix.(control_M, 0.0; force=true)
    JuMP.fix.(control_F[:, :, wildtype, :, :], 0.0; force=true)
    JuMP.fix.(control_F[:, :, :, homozygous_modified, :], 0.0; force=true)

    return
end

struct _MixedStrategy{T, U} end

function _add_constraint(
    model::JuMP.Model,
    do_binary::Bool,
    optinfo_dict::Dict,
    node_number::Int64,
    life_stage::Tuple,
    node_strategy::ReleaseStrategy,
    homozygous_modified,
    wildtype,
)
    for entry in life_stage
        if !(entry <: LifeStage)
            error("$entry is not a LifeStage.")
        end
    end

    life_stage = _MixedStrategy{life_stage[1], life_stage[2]}

    if life_stage == _MixedStrategy{Male, Female} ||
       life_stage == _MixedStrategy{Female, Male}
        _add_constraint(
            model,
            do_binary::Bool,
            optinfo_dict,
            node_number,
            _MixedStrategy{Male, Female},
            node_strategy,
            homozygous_modified,
            wildtype,
        )
    else
        error("The release strategy $(life_stage) is not functional.")
    end
end

function _add_constraint(
    model::JuMP.Model,
    do_binary::Bool,
    optinfo_dict::Dict,
    node_number::Int64,
    life_stage::Type{_MixedStrategy{Male, Female}},
    node_strategy::ReleaseStrategy,
    homozygous_modified,
    wildtype,
)
    data = optinfo_dict["constraints"]
    max_per_timestep = node_strategy.release_size_max_per_timestep
    min_per_timestep = node_strategy.release_size_min_per_timestep
    max_overall = node_strategy.release_max_over_timehorizon #TODO: fix
    control_F = model[:control_F]
    control_M = model[:control_M]
    release_location = model[:release_location]
    sets = model.obj_dict[:Sets]
    N = sets[:N]
    O = sets[:O]
    SF = sets[:SF]
    SM = sets[:SM]
    G = sets[:G]
    T = sets[:T]

    CONTAINER_control_limit_schedule_upper_F =
        model.obj_dict[:control_limit_schedule_upper_F]
    CONTAINER_control_limit_schedule_lower_F =
        model.obj_dict[:control_limit_schedule_lower_F]
    CONTAINER_control_limit_total_F = model.obj_dict[:control_limit_total_F]

    CONTAINER_control_limit_schedule_upper_M =
        model.obj_dict[:control_limit_schedule_upper_M]
    CONTAINER_control_limit_schedule_lower_M =
        model.obj_dict[:control_limit_schedule_lower_M]
    CONTAINER_control_limit_total_M = model.obj_dict[:control_limit_total_M]

    CONTAINER_control_equivalence = model.obj_dict[:control_equivalence]

    range =
        (node_strategy.release_start_time):(node_strategy.release_time_interval):(node_strategy.release_end_time)

    # EQUIVALENCE
    for n in node_number, o in O, t in T
        CONTAINER_control_equivalence[n, o, t] = JuMP.@constraint(
            model,
            control_F[n, o, homozygous_modified, wildtype, t] ==
            control_M[n, o, 1, homozygous_modified, t]
        )
    end

    for n in node_number, t in T
        if t ∈ range
            if do_binary
                #@info "made binary for node $n and time $t"
                JuMP.set_binary(release_location[n, t])
            end

            for o in O, g in G
                max_v = g == homozygous_modified ? max_per_timestep : 0.0
                min_v = g == homozygous_modified ? min_per_timestep : 0.0

                # SCHEDULE F
                for s in SF
                    #@info "max used in $s, $g is FEMALES = $max_v"
                    #@info "min used in $s, $g is FEMALES = $min_v"

                    CONTAINER_control_limit_schedule_upper_F[n, o, s, g, t] =
                        JuMP.@constraint(
                            model,
                            # max
                            control_F[n, o, g, s, t] <= release_location[n, t] * max_v
                        )
                    CONTAINER_control_limit_schedule_lower_F[n, o, s, g, t] =
                        JuMP.@constraint(
                            model,
                            # min
                            control_F[n, o, g, s, t] >= release_location[n, t] * min_v
                        )
                end

                # SCHEDULE M
                for s in SM
                    #@info "max used in $s, $g is MALES = $max_v"
                    #@info "min used in $s, $g is MALES = $min_v"

                    CONTAINER_control_limit_schedule_upper_M[n, o, s, g, t] =
                        JuMP.@constraint(
                            model,
                            # max
                            control_M[n, o, s, g, t] <= release_location[n, t] * max_v
                        )
                    CONTAINER_control_limit_schedule_lower_M[n, o, s, g, t] =
                        JuMP.@constraint(
                            model,
                            # min
                            control_M[n, o, s, g, t] >= release_location[n, t] * min_v
                        )
                end
            end
        else
            JuMP.fix(release_location[n, t], 0.0; force=true)
            for o in O, g in G
                for s in SF
                    #@info "fixed FEMALES to 0.0 in node $n, time $t"
                    CONTAINER_control_limit_schedule_upper_F[n, o, s, g, t] =
                        JuMP.@constraint(model, control_F[n, o, g, s, t] <= 0)
                    CONTAINER_control_limit_schedule_lower_F[n, o, s, g, t] =
                        JuMP.@constraint(model, control_F[n, o, g, s, t] >= 0)
                end
                for s in SM
                    #@info "fixed MALES to 0.0 in node $n, time $t"
                    CONTAINER_control_limit_schedule_upper_M[n, o, s, g, t] =
                        JuMP.@constraint(model, control_M[n, o, s, g, t] <= 0)
                    CONTAINER_control_limit_schedule_lower_M[n, o, s, g, t] =
                        JuMP.@constraint(model, control_M[n, o, s, g, t] >= 0)
                end
            end
        end
    end

    # MAXIMUM F
    @info "overall max FEMALES = $(max_overall*MALE_FEMALE_RELEASE_FRACTION)"
    CONTAINER_control_limit_total_F = JuMP.@constraint(
        model,
        sum(control_F) <= max_overall * MALE_FEMALE_RELEASE_FRACTION
    )

    # MAXIMUM M
    @info "overall max MALES = $(max_overall*MALE_FEMALE_RELEASE_FRACTION)"
    CONTAINER_control_limit_total_M = JuMP.@constraint(
        model,
        sum(control_M) <= max_overall * MALE_FEMALE_RELEASE_FRACTION
    )

    JuMP.fix.(control_M[:, :, :, wildtype, :], 0.0; force=true)
    JuMP.fix.(control_F[:, :, wildtype, :, :], 0.0; force=true)
    JuMP.fix.(control_F[:, :, :, homozygous_modified, :], 0.0; force=true)

    return
end

function _populate_constraints_dict(
    releaseconstraints_dict,
    node_number,
    strategy,
    life_stage::Type{<:LifeStage},
)
    for tx in
        (strategy.release_start_time):(strategy.release_time_interval):(strategy.release_end_time)
        releaseconstraints_dict[node_number]["org_con"][1]["stage_con"][life_stage][strategy.release_this_gene_index][tx] =
            strategy.release_size_max_per_timestep
    end
    return
end

function _populate_constraints_dict(
    releaseconstraints_dict,
    node_number,
    strategy,
    life_stage::Tuple,
)
    for entry in life_stage
        if !(entry <: LifeStage)
            error("$entry is not a LifeStage.")
        end
    end

    lifestage_entry = _MixedStrategy{life_stage[1], life_stage[2]}
    _populate_constraints_dict(
        releaseconstraints_dict,
        node_number,
        strategy,
        lifestage_entry,
    )
end

function _populate_constraints_dict(
    releaseconstraints_dict,
    node_number,
    strategy,
    ::Type{_MixedStrategy{Male, Female}},
)
    range =
        (strategy.release_start_time):(strategy.release_time_interval):(strategy.release_end_time)
    for tx in range
        releaseconstraints_dict[node_number]["org_con"][1]["stage_con"][Male][strategy.release_this_gene_index][tx] =
            strategy.release_size_max_per_timestep * MALE_FEMALE_RELEASE_FRACTION
        releaseconstraints_dict[node_number]["org_con"][1]["stage_con"][Female][strategy.release_this_gene_index][tx] =
            strategy.release_size_max_per_timestep * MALE_FEMALE_RELEASE_FRACTION
    end
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

    # Create empty constraints dict
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
            @warn "Release end time specified in node $(node_number) occurs outside the time horizon (tspan = $(tspan))."
        end

        if strategy.release_size_min_per_timestep >= strategy.release_size_max_per_timestep
            @warn "Release size minimum per timestep specified in node $(node_number) is greater than or equal to release size maximum per timestep."
        end

        if all(
            isequal(
                node_strategy[node_number].release_max_over_timehorizon,
                strategy.release_max_over_timehorizon,
            ),
        ) == false
            @warn "Entries for the `ReleaseStrategy` field `release_max_over_timehorizon` do not match; this field should be equivalent for all nodes in the network."
        end

        strategy.release_start_time === nothing ? (strategy.release_start_time = tspan[1]) :
        strategy.release_start_time
        strategy.release_end_time === nothing ? (strategy.release_end_time = tspan[2]) :
        strategy.release_end_time

        # Fill release dict
        ####################
        _populate_constraints_dict(
            releaseconstraints_dict,
            node_number,
            node_strategy[node_number],
            node_strategy[node_number].release_this_life_stage,
        )
    end

    # Include release dict
    ####################
    optinfo_dict["constraints"] = releaseconstraints_dict

    # Add constraints to model
    ####################
    for (node_number, strategy) in enumerate(node_strategy)
        _add_constraint(
            model,
            do_binary::Bool,
            optinfo_dict,
            node_number,
            node_strategy[node_number].release_this_life_stage,
            node_strategy[node_number],
            homozygous_modified,
            wildtype,
        )
    end

    return optinfo_dict
end
