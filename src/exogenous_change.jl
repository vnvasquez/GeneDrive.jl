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
    function ExogenousInputs(; # ; for default arguments
        intervention = Dict{Symbol, Dict}(),
        temperature = Dict{Symbol, Float64}())
        new(intervention,
        temperature)
    end
end

function _make_intervention_space(::Type{<:LifeStage}, number_of_genes::Int)
    return zeros(Float64, number_of_genes)
end

function _make_intervention_space(::Type{Female}, number_of_genes::Int)
    return zeros(Float64, number_of_genes, number_of_genes)
end

function _allocate_intervention_space(node::Node) # _ for internal function
        #upper level dict
        values = Dict{Symbol, Dict}(get_name(node) => Dict())
        # second level dict
        for organism in get_organisms(node)
            values[get_name(node)][organism] = Dict()
            number_of_genes = count_genotypes(node, organism)
            for life_stage in keys(get_lifestages(node, organism))
                #third level dict
                values[get_name(node)][organism][life_stage] = _make_intervention_space(life_stage, number_of_genes)
            end
        end
        return values
end

function _get_temperature_value(node::Node)
   return Dict(get_name(node) => initialize_temperature_model(node.temperature)) # second argument = initial time @ 0.0 returns initial temp
end

"""
        function ExogenousInputs(node::Node)

    Returns allocated space for exogenous inputs to a `Node`.
"""
function ExogenousInputs(node::Node)
    return ExogenousInputs(intervention = _allocate_intervention_space(node),
    temperature = _get_temperature_value(node))
end

"""
        function ExogenousInputs(network::Network)

    Returns allocated space for exogenous inputs to a `Network`.
"""
function ExogenousInputs(network::Network)
    values_intervention = Dict{Symbol, Dict}()
    values_temp = Dict{Symbol, Float64}()
    for (name, node) in get_nodes(network)
        merge!(values_intervention, _allocate_intervention_space(node))
        merge!(values_temp, _get_temperature_value(node))
    end
    return ExogenousInputs(intervention = values_intervention, temperature = values_temp)
end

"""
        function get_exogenous_intervention(inputs::ExogenousInputs, node::Node, organism, stage, gene)

    Returns biological control intervention relevant to specified `Node`, `Organism`, `LifeStage`, and `Genotype`.
"""
function get_exogenous_intervention(inputs::ExogenousInputs, node::Node, organism, stage, gene)
    node_name = get_name(node)
    return inputs.intervention[node_name][organism][stage][gene]
end

"""
        function get_exogenous_intervention(inputs::ExogenousInputs, node::Node, organism, stage, gene)

    Returns biological control intervention relevant to specified `Female` life stage.
"""
function get_exogenous_intervention(inputs::ExogenousInputs, node::Node, organism, stage::Type{Female}, gene, female_index)
    node_name = get_name(node)
    return inputs.intervention[node_name][organism][stage][female_index, gene]
end

"""
        function get_exogenous_temperature(inputs::ExogenousInputs, node)

    Returns temperature relevant to specified `Node`.
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
            set::diffeq.CallbackSet
        end

    Data for biological control interventions predicated on releasing modified organisms. Applicable to both suppression and replacement techniques.

# Arguments
- `node::Symbol`: `Node` where releases occur.
- `organism::Type{<:Species}`: Species of organism being released.
- `stage::Type{<:LifeStage}`: `Life stage` of organism being released.
- `gene_index::Int64`: Genotype being released to implement the intervention.
- `times::Vector{Float64}`: Release start/stop time points. Multiple sets of start/stop points with variable intervals permitted.
- `values::Vector{Float64}`: Number of modified organisms released during each time period. Variably sized releases permitted.
- `set::diffeq.CallbackSet`
"""
mutable struct Release <: Intervention
    node::Symbol
    organism::Type{<:Species}
    stage::Type{<:LifeStage}
    gene_index::Int64 #TODO: consider updating to string or type rather than number.
    times::Vector{Float64}
    values::Vector{Int64}
    set::diffeq.CallbackSet
end

"""
        function Release(node::Node, organism, stage::Type{T}, gene, times::Vector, variable_release::Vector{Float64}) where T <: LifeStage

    Returns `Release` object specifying details of biological control intervention where release size is variable over time.
"""
function Release(node::Node, organism, stage::Type{T}, gene, times::Vector, variable_release::Vector{Int64}) where T <: LifeStage
    @assert length(times) == length(variable_release)
    callbacks = Vector()
    stop_times = Vector()
    for (index, time) in enumerate(times)
        condition_start = (u, t, integrator) -> t == time
        affect_start = get_effect_intervention(node, organism, stage, gene, variable_release[index])
        push!(callbacks, diffeq.DiscreteCallback(condition_start, affect_start))
        stop_time = time + 1.01 #TODO: Improve this stop_times approach: 1.01*integrator.dt
        condition_stop = (u, t, integrator) -> t == stop_time #time + integrator.dt
        push!(stop_times, stop_time)
        affect_stop = get_effect_intervention(node, organism, stage, gene, 0.0)
        push!(callbacks, diffeq.DiscreteCallback(condition_stop, affect_stop))
    end
    callback_set = diffeq.CallbackSet((), tuple(callbacks...))
    return Release(get_name(node), organism, stage, gene, vcat(times, stop_times), variable_release, callback_set)
end

"""
        function Release(node::Node, organism, stage, gene, times::Vector, fixed_release::Float64)

    Returns `Release` object specifying details of biological control intervention where release size is fixed over time.
"""
function Release(node::Node, organism, stage::Type{T}, gene, times::Vector, fixed_release::Int64) where T <: LifeStage
    release_values = fill(fixed_release, length(times))
    return Release(node, organism, stage, gene, times, release_values)
end

function get_effect_intervention(node::Node, organism, stage::Type{<:LifeStage}, gene, value)
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
        function TemperatureSeriesData(node::Node, times::Vector, values::Vector)

    Returns `TemperatureSeriesData` instance populated with specified data for use in simulation.
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
    times::Vector{Tuple{Float64,Float64}}
    values::Vector{Float64}
    set::diffeq.CallbackSet
end

"""
        function TemperatureShockData(node::Node, times::Vector{Tuple{Float64,Float64}}, variable_shock::Vector{Float64})

    Returns `TemperatureShockData` instance specifying data for variable shock values.
"""
function TemperatureShockData(node::Node, times::Vector{Tuple{Float64,Float64}}, variable_shock::Vector{Float64})
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

    Returns `TemperatureShockData` instance specifying data for fixed (all the same) shock values.
"""
function TemperatureShockData(node::Node, times::Vector{Tuple{Float64,Float64}}, fixed_shock::Float64)
    shock_values = fill(fixed_shock, length(times))
    return TemperatureShockData(node, times, shock_values)
end

function get_effect_temperature(node_name::Symbol, value::Float64)
    return (integrator) -> begin
        integrator.p[2].temperature[node_name] = value
    end
end
