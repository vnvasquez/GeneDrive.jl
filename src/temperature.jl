################################################################################
#                           Climatic Temperature                               #
################################################################################

abstract type Temperature end

########################################
#          ConstantTemperature         #
########################################
"""
    mutable struct ConstantTemperature <: Temperature end

Data for simulation that features a single constant temperature in °C. Generated internally.

# Arguments

  - `value::Float64`: Temperature value in °C; re-used by model at every time step.
"""
mutable struct ConstantTemperature <: Temperature
    value::Float64
end

"""
    get_temperature_value(temperature_model::ConstantTemperature, temp_value_from_inputs::Float64, t)

Return value of temperature in °C for specifed time step. Accounts for perturbations to baseline temperature model where applicable.
"""
function get_temperature_value(
    temperature_model::ConstantTemperature,
    temp_value_from_inputs::Float64,
    t,
)
    return temperature_model.value + temp_value_from_inputs
end

"""
    initialize_temperature_model(data::ConstantTemperature)

Return value of temperature in °C in first time step, for simulation initialization.
"""
function initialize_temperature_model(data::ConstantTemperature)
    # time[0] = time[:]
    return data.value
end

########################################
#        SinusoidalTemperature         #
########################################
"""
    mutable struct SinusoidalTemperature <: Temperature
        a::Float64
        b::Int64
        c::Int64
        d::Float64
    end

Data for simulation that features sinusoidal temperature fluctuation in °C. Uses `cosine` implementation applicable to the Southern Hemisphere. Generated internally.

# Arguments

  - `a::Float64`: Amplitude.
  - `b::Float64`: Periodicity coefficient.
  - `c::Float64`: Time period (days).
  - `d::Float64`: Mean.
"""
mutable struct SinusoidalTemperature <: Temperature
    a::Float64
    b::Int64
    c::Int64
    d::Float64
end

"""
    get_temperature_value(temperature_model::SinusoidalTemperature, temp_value_from_inputs::Float64, t)

Return value of temperature in °C for specifed time step. Accounts for perturbations to baseline temperature model where applicable.
"""
function get_temperature_value(
    temperature_model::SinusoidalTemperature,
    temp_value_from_inputs::Float64,
    t,
)
    # NB: temperature_model.d sets initial value;
    # inputs used thereafter to account for mean because allows for temperature_model.d + shocks as needed
    return temperature_model.a * cos((temperature_model.b * π / temperature_model.c) * t) +
           temp_value_from_inputs
end

"""
    initialize_temperature_model(data::SinusoidalTemperature)

Return value of temperature in °C in first time step, for simulation initialization.
"""
function initialize_temperature_model(data::SinusoidalTemperature)
    # time = 0.0
    return data.a * cos((data.b * π / data.c) * 0.0) + data.d
end

########################################
#        TimeSeriesTemperature         #
########################################
"""
    mutable struct TimeSeriesTemperature <: Temperature end

Data for simulation that uses temperature time series in °C.

# Arguments

  - `value::Float64`: Time series of temperature values in °C.
"""
mutable struct TimeSeriesTemperature <: Temperature
    values::Vector{Float64}
end

"""
    get_temperature_value(temperature_model::TimeSeriesTemperature, temp_value_from_inputs::Float64, t)

Return value of temperature in °C for specifed time step. Perturbations to baseline temperature, where applicable, should be directly added to time series prior to running model.
"""
function get_temperature_value(
    temperature_model::TimeSeriesTemperature,
    temp_value_from_inputs::Float64,
    t,
)
    return temp_value_from_inputs
end

"""
    initialize_temperature_model(data::TimeSeriesTemperature)

Return first value of temperature in °C. Used for simulation initialization.
"""
function initialize_temperature_model(data::TimeSeriesTemperature)
    # value[1] = time[0]
    return data.values[1]
end

########################################
#         ScenarioTemperature          #
########################################
"""
    mutable struct ScenarioTemperature <: Temperature end

Data for simulation that uses ensembles of temperature time series in °C.

# Arguments

  - `data::Matrix{Float64}`: Array of ensemble members, each column of which is a time series of temperature values in °C.
  - `probability::Vector{Float64}`: Vector of probabilities with which given scenarios are expected to occur (one probability applicable per ensemble member).
  - `selected_scenario::Int`: Index to run dynamic model over selected scenarios; `Int` argument required. For decision model, use `nothing`.
"""
mutable struct ScenarioTemperature <: Temperature
    values::Matrix{Float64}
    probability::Vector{Float64}
    selected_scenario::Union{Nothing, Int}

   # function ScenarioTemperature(
   #     values::Matrix{Float64},
   #     probability::Vector{Float64})
   #    new(values, probability, nothing)
   # end
end

function ScenarioTemperature(
    values::Matrix{Float64},
    probability::Vector{Float64};
    selected_scenario=nothing)
    return ScenarioTemperature(values, probability, selected_scenario)
end

"""
    get_temperature_value(temperature_model::ScenarioTemperature, temp_value_from_inputs::Float64, t)

Return value of temperature in °C for selected scenario and specifed time step. Perturbations to baseline temperature, where applicable, should be directly added to time series prior to running model.
"""
function get_temperature_value(
    temperature_model::ScenarioTemperature,
    temp_value_from_inputs::Float64,
    t,
)
    return temp_value_from_inputs
end

"""
    initialize_temperature_model(data::ScenarioTemperature)

Return first value of temperature in °C for selected scenario. Used for simulation initialization.
"""
function initialize_temperature_model(data::ScenarioTemperature)
    return data.values[data.selected_scenario][1]
end
