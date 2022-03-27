################################################################################
#                                    General                                   #
################################################################################

"""
    struct ExogenousInputs
        intervention::Dict{Symbol, Dict}
        temperature::Dict{Symbol, Float64}
        function ExogenousInputs(;
            intervention = Dict{Symbol, Dict}(),
            temperature = Dict{Symbol, Float64}())
            new(intervention,
            temperature)
        end
    end

Data for perturbations to the simulated system, including (1) biological control interventions and (2) externally defined temperature series/shocks.
"""
struct ExogenousInputs
    intervention::Dict{Symbol, Dict}
    temperature::Dict{Symbol, Float64}
    function ExogenousInputs(;
        intervention=Dict{Symbol, Dict}(),
        temperature=Dict{Symbol, Float64}(),
    )
        new(intervention, temperature)
    end
end

function _make_intervention_space(::Type{<:LifeStage}, number_of_genes::Int)
    return zeros(Float64, number_of_genes)
end

function _make_intervention_space(::Type{Female}, number_of_genes::Int)
    return zeros(Float64, number_of_genes, number_of_genes)
end

function _allocate_intervention_space(node::Node)
    values = Dict{Symbol, Dict}(get_name(node) => Dict())
    for organism in get_organisms(node)
        values[get_name(node)][organism] = Dict()
        number_of_genes = count_genotypes(node, organism)
        for life_stage in keys(get_lifestages(node, organism))
            values[get_name(node)][organism][life_stage] =
                _make_intervention_space(life_stage, number_of_genes)
        end
    end
    return values
end

function _get_temperature_value(node::Node)
    return Dict(get_name(node) => initialize_temperature_model(node.temperature)) # second argument = initial time @ 0.0 returns initial temp
end

"""
    ExogenousInputs(node::Node)

Return allocated space for exogenous inputs to `Node`.
"""
function ExogenousInputs(node::Node)
    return ExogenousInputs(
        intervention=_allocate_intervention_space(node),
        temperature=_get_temperature_value(node),
    )
end

"""
    ExogenousInputs(network::Network)

Return allocated space for exogenous inputs to a `Network`.
"""
function ExogenousInputs(network::Network)
    values_intervention = Dict{Symbol, Dict}()
    values_temp = Dict{Symbol, Float64}()
    for (name, node) in get_nodes(network)
        merge!(values_intervention, _allocate_intervention_space(node))
        merge!(values_temp, _get_temperature_value(node))
    end
    return ExogenousInputs(intervention=values_intervention, temperature=values_temp)
end

"""
        function get_exogenous_intervention(inputs::ExogenousInputs, node::Node, organism, stage, gene)

    Returns biological control intervention relevant to specified `Node`, `Organism`, `LifeStage`, and `Genotype`.
"""
function get_exogenous_intervention(
    inputs::ExogenousInputs,
    node::Node,
    organism,
    stage,
    gene,
)
    node_name = get_name(node)
    return inputs.intervention[node_name][organism][stage][gene]
end

"""
    get_exogenous_intervention(inputs::ExogenousInputs, node::Node, organism, stage::Type{Female}, gene, female_index)

Return biological control intervention relevant to `Female`.
"""
function get_exogenous_intervention(
    inputs::ExogenousInputs,
    node::Node,
    organism,
    stage::Type{Female},
    gene,
    female_index,
)
    node_name = get_name(node)
    return inputs.intervention[node_name][organism][stage][female_index, gene]
end

"""
    get_exogenous_temperature(inputs::ExogenousInputs, node)

Return temperature relevant to `Node`.
"""
function get_exogenous_temperature(inputs::ExogenousInputs, node)
    return inputs.temperature[node]
end

abstract type ExogenousChange end

################################################################################
#                                  Interventions                               #
################################################################################

abstract type Intervention <: ExogenousChange end

"""
    mutable struct Release <: Intervention
        node::Symbol
        organism::Type{<:Species}
        stage::Type{<:LifeStage}
        gene_index::Int64
        times::Vector{Float64}
        values::Vector{Float64}
        callbacks::Vector
    end

Data for biological control interventions predicated on releasing modified organisms. Applicable to both suppression and replacement techniques.

# Arguments

  - `node::Symbol`: `Node` where releases occur.
  - `organism::Type{<:Species}`: Species of organism being released.
  - `stage::Type{<:LifeStage}`: `Life stage` of organism being released.
  - `gene_index::Int64`: Genotype being released to implement the intervention.
  - `times::Vector{Float64}`: Release start/stop time points. Multiple sets of start/stop points with variable intervals permitted.
  - `values::Vector{Float64}`: Number of modified organisms released during each time period. Variably sized releases permitted.
  - `callbacks::Vector`
"""
mutable struct Release <: Intervention
    node::Symbol
    organism::Type{<:Species}
    stage::Type{<:LifeStage}
    gene_index::Int64
    times::Vector{Float64}
    values::Vector{Int64}
    callbacks::Vector
end

"""
    Release(node::Node, organism, stage::Type{T}, gene, times::Vector, variable_release::Vector{Float64}) where T <: LifeStage

Return `Release` object specifying details of biological control intervention where release size is variable over time.
"""
function Release(
    node::Node,
    organism,
    stage::Type{T},
    gene,
    times::Vector,
    variable_release::Vector{Int64},
) where {T <: LifeStage}
    @assert length(times) == length(variable_release)
    callbacks = Vector()
    stop_times = Vector()
    for (index, time) in enumerate(times)
        condition_start = (u, t, integrator) -> t == time
        affect_start =
            get_effect_intervention(node, organism, stage, gene, variable_release[index])
        push!(callbacks, diffeq.DiscreteCallback(condition_start, affect_start))
        stop_time = time + 1.01 #TODO: Improve this stop_times approach: 1.01*integrator.dt
        condition_stop = (u, t, integrator) -> t == stop_time #time + integrator.dt
        push!(stop_times, stop_time)
        affect_stop = get_effect_intervention(node, organism, stage, gene, 0.0)
        push!(callbacks, diffeq.DiscreteCallback(condition_stop, affect_stop))
    end
    return Release(
        get_name(node),
        organism,
        stage,
        gene,
        vcat(times, stop_times),
        variable_release,
        callbacks,
    )
end

"""
    Release(node::Node, organism, stage, gene, times::Vector, fixed_release::Float64)

Return `Release` object specifying details of biological control intervention where release size is fixed over time.
"""
function Release(
    node::Node,
    organism,
    stage::Type{T},
    gene,
    times::Vector,
    fixed_release::Int64,
) where {T <: LifeStage}
    release_values = fill(fixed_release, length(times))
    return Release(node, organism, stage, gene, times, release_values)
end

function get_effect_intervention(
    node::Node,
    organism,
    stage::Type{<:LifeStage},
    gene,
    value,
)
    node_name = get_name(node)
    return (integrator) -> begin
        integrator.p[2].intervention[node_name][organism][stage][gene] = value
    end
end

function get_effect_intervention(node::Node, organism, stage::Type{Female}, gene, value)
    wild_type = get_wildtype(node, organism)
    node_name = get_name(node)
    return (integrator) -> begin
        integrator.p[2].intervention[node_name][organism][stage][gene, wild_type] = value
    end
end

"""
    mutable struct ProportionalRelease <: Intervention
        node::Symbol
        organism::Type{<:Species}
        stage::Type{<:LifeStage}
        gene_index::Int64
        times::Vector{Float64}
        values::Vector{Float64}
        callbacks::Vector
        adults_counting::String
    end

Data for biological control interventions predicated on releasing modified organisms. Applicable to both suppression and replacement techniques.

# Arguments

  - `node::Symbol`: `Node` where releases occur.
  - `organism::Type{<:Species}`: Species of organism being released.
  - `stage::Type{<:LifeStage}`: `Life stage` of organism being released.
  - `gene_index::Int64`: Genotype being released to implement the intervention.
  - `times::Vector{Float64}`: Release start/stop time points. Multiple sets of start/stop points with variable intervals permitted.
  - `values::Vector{Float64}`: Number of modified organisms released during each time period. Variably sized releases permitted.
  - `callbacks::Vector`
  - `adults_to_count::String`: The population against which to measure proportional release size. Choices include `Male`, `Female`, or `All`.
"""
mutable struct ProportionalRelease <: Intervention
    node::Symbol
    organism::Type{<:Species}
    stage::Type{<:LifeStage}
    gene_index::Int64
    times::Vector{Float64}
    proportion::Float64
    callbacks::Vector
    adults_to_count::String
end

function _get_adult_map(counting_string)
    if counting_string == "Females"
        return (1.0, 0.0)
    elseif counting_string == "Males"
        return (0.0, 1.0)
    elseif counting_string == "All"
        return (1.0, 1.0)
    else
        error(
            "$counting_string not recognized. Only `Females`, `Males`, or `All` are accepted.",
        )
    end
end

"""
    ProportionalRelease(node::Node, organism, stage::Type{T}, gene, times::Vector, release::Float64, adult_counting::String) where T <: LifeStage

Return `ProportionalRelease` object specifying details of biological control intervention where release size is predicated on real-time wild population level.
"""
function ProportionalRelease(
    node::Node,
    organism,
    stage::Type{T},
    gene,
    times::Vector,
    release::Float64,
    adult_counting::String,
) where {T <: LifeStage}
    callbacks = Vector()
    stop_times = Vector()
    for (index, time) in enumerate(times)
        condition_start = (u, t, integrator) -> t == time
        affect_start = get_effect_intervention_proportional(
            node,
            organism,
            stage,
            gene,
            release,
            _get_adult_map(adult_counting),
        )
        push!(callbacks, diffeq.DiscreteCallback(condition_start, affect_start))
        stop_time = time + 1.01 #TODO: Improve this stop_times approach: 1.01*integrator.dt
        condition_stop = (u, t, integrator) -> t == stop_time #time + integrator.dt
        push!(stop_times, stop_time)
        affect_stop = get_effect_intervention(node, organism, stage, gene, 0.0)
        push!(callbacks, diffeq.DiscreteCallback(condition_stop, affect_stop))
    end
    return ProportionalRelease(
        get_name(node),
        organism,
        stage,
        gene,
        vcat(times, stop_times),
        release,
        callbacks,
        adult_counting,
    )
end

function release_value(current_size, value)
    return clamp(current_size * value, 10, 100000)
end

function get_effect_intervention_proportional(
    node::Node,
    organism,
    stage::Type{Male},
    gene,
    value,
    adults_map,
)
    wild_type = get_wildtype(node, organism)
    node_name = get_name(node)
    number_of_genes = length(get_genotypes(node, organism))
    index_males =
        sum(count_substages(node, organism)[1:(end - 1)]):(sum(
            count_substages(node, organism)[1:end],
        ) - 1)
    index_just_females = sum(count_substages(node, organism)) + number_of_genes - 1
    return (integrator) -> begin
        # Hardcoding the [1] in the x don't have node indexing yet (TODO: come back to this when focused on network)
        current_size_males = sum(integrator.u.x[1][index_males, wild_type])
        current_size_females = sum(integrator.u.x[1][index_just_females, :])
        current_size =
            adults_map[1] * current_size_females + adults_map[2] * current_size_males
        release_value(current_size, value)
        integrator.p[2].intervention[node_name][organism][stage][gene] =
            release_value(current_size, value)
    end
end

function get_effect_intervention_proportional(
    node::Node,
    organism,
    stage::Type{Female},
    gene,
    value,
    adults_map,
)
    wild_type = get_wildtype(node, organism)
    node_name = get_name(node)
    number_of_genes = length(get_genotypes(node, organism))
    index_males =
        sum(count_substages(node, organism)[1:(end - 1)]):(sum(
            count_substages(node, organism)[1:end],
        ) - 1)
    index_just_females = sum(count_substages(node, organism)) + number_of_genes - 1
    return (integrator) -> begin
        # Hardcoding the [1] in the x don't have node indexing yet (TODO: come back to this when focused on network)
        current_size_males = sum(integrator.u.x[1][index_males, wild_type])
        current_size_females = integrator.u.x[1][index_just_females, wild_type]
        current_size =
            adults_map[1] * current_size_females + adults_map[2] * current_size_males
        release_value(current_size, value)
        integrator.p[2].intervention[node_name][organism][stage][gene, wild_type] =
            release_value(current_size, value)
    end
end

################################################################################
#                                  Temperature                                 #
################################################################################

"""
    mutable struct TemperatureSeriesData <: ExogenousChange
        node::Symbol
        times::Vector{Float64}
        values::Vector{Float64}
        set::diffeq.CallbackSet
    end

Data for time series of temperature.

# Arguments

  - `node::Symbol`: `Node` where temperature is realized.
  - `times::Vector{Float64}`: Time step (day) of realized temperature value.
  - `values::Vector{Float64}`: Temperature in °C.
  - `set::diffeq.CallbackSet`
"""
mutable struct TemperatureSeriesData <: ExogenousChange
    node::Symbol
    times::Vector{Float64}
    values::Vector{Float64}
    set::diffeq.CallbackSet
end

"""
        TemperatureSeriesData(node::Node, times::Vector, values::Vector)

Return `TemperatureSeriesData` instance populated with specified data for use in simulation.
"""

function TemperatureSeriesData(node::Node, times::Vector, values::Vector)
    @assert length(times) == length(values)
    callbacks = Vector()
    for (index, time) in enumerate(times)
        condition = (u, t, integrator) -> t == time
        affect = get_effect_temperature(get_name(node), values[index])
        push!(callbacks, diffeq.DiscreteCallback(condition, affect))
    end
    callback_set = diffeq.CallbackSet((), tuple(callbacks...))
    return TemperatureSeriesData(get_name(node), times, values, callback_set)
end

"""
    mutable struct TemperatureShockData <: ExogenousChange
        node::Symbol
        times::Vector{Float64}
        values::Vector{Float64}
        set::diffeq.CallbackSet
    end

Data for time series of temperature.

# Arguments

  - `node::Symbol`: `Node` where temperature shock occurs.
  - `times::Vector{Tuple{Float64,Float64}}`: Vector of tuples defining shock start/stop time points. Multiple discrete shock periods may be included.
  - `values::Vector{Float64}`: Single temperature value in °C added to the baseline temperature during the specified time period. Multiple discrete temperature changes may be included; values may be positive or negative.
  - `set::diffeq.CallbackSet`
"""
mutable struct TemperatureShockData <: ExogenousChange
    node::Symbol
    times::Vector{Tuple{Float64, Float64}}
    values::Vector{Float64}
    set::diffeq.CallbackSet
end

"""
    TemperatureShockData(node::Node, times::Vector{Tuple{Float64,Float64}}, variable_shock::Vector{Float64})

Return `TemperatureShockData` instance specifying data for variable shock values.
"""
function TemperatureShockData(
    node::Node,
    times::Vector{Tuple{Float64, Float64}},
    variable_shock::Vector{Float64},
)
    @assert length(times) == length(variable_shock)
    callbacks = Vector()
    for (index, time) in enumerate(times)
        condition_start = (u, t, integrator) -> t == time[1]
        affect_start = get_effect_temperature(get_name(node), variable_shock[index])
        push!(callbacks, diffeq.DiscreteCallback(condition_start, affect_start))
        condition_stop = (u, t, integrator) -> t == time[2]
        affect_stop = get_effect_temperature(get_name(node), 0.0)
        push!(callbacks, diffeq.DiscreteCallback(condition_stop, affect_stop))
    end
    callback_set = diffeq.CallbackSet((), tuple(callbacks...))
    times = collect(Iterators.flatten(times))
    return TemperatureShockData(get_name(node), times, variable_shock, callback_set)
end

"""
    function TemperatureShockData(node::Node, times::Vector{Tuple{Float64,Float64}}, fixed_shock::Float64)

Return `TemperatureShockData` instance specifying data for fixed (all the same) shock values.
"""
function TemperatureShockData(
    node::Node,
    times::Vector{Tuple{Float64, Float64}},
    fixed_shock::Float64,
)
    shock_values = fill(fixed_shock, length(times))
    return TemperatureShockData(node, times, shock_values)
end

function get_effect_temperature(node_name::Symbol, value::Float64)
    return (integrator) -> begin
        integrator.p[2].temperature[node_name] = value
    end
end
