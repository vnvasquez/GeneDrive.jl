
################################################################################
#                              Density Dependence                              #
################################################################################

abstract type DensityDependence end

"""
    mutable struct Density{D <: DensityDependence}
        model::Type{D}
        param::Float64
    end

Data for density dependence model.

# Arguments

  - `model::Type{D}`: User-selected density dependence formulation. Includes `LogisticDensity`, `LinearDensity`, and `NoDensity`.
  - `param::Float64`: Calculated internally by `init` function.
"""
mutable struct Density{D <: DensityDependence}
    model::Type{D}
    param::Float64
end

"""
    struct LogisticDensity <: DensityDependence end

Density dependence model of type `LogisticDensity`.
"""
struct LogisticDensity <: DensityDependence end

"""
    compute_density(data::Density{LogisticDensity}, stage)

Return the effect of `LogisticDensity` on life stage instance. Calculated each time step.
"""
function compute_density(data::Density{LogisticDensity}, stage)
    return (1 + (sum(stage) / data.param))
end

"""
    struct LinearDensity <: DensityDependence end

Density dependence model of type `LinearDensity`.
"""
struct LinearDensity <: DensityDependence end

"""
    compute_density(data::Density{LinearDensity}, stage)

Return the effect of `LinearDensity` on life stage instance. Calculated each time step.
"""
function compute_density(data::Density{LinearDensity}, stage)
    return data.param * sum(stage)
    # TODO: check against Hancock paper.
end

"""
    struct NoDensity <: DensityDependence end

Density dependence model of type `NoDensity`.
"""
struct NoDensity <: DensityDependence end

"""
    compute_density(data::Density{NoDensity}, stage)

Return the effect of `NoDensity` on life stage instance. Calculated each time step.
"""
function compute_density(data::Density{NoDensity}, stage)
    return data.param
end

################################################################################
#                                Life Stages                                   #
################################################################################

abstract type LifeStage end

"""
    mutable struct Stage{L <: LifeStage}
        μ_temperature_response::TemperatureResponse
        q_temperature_response::TemperatureResponse
        n::Union{Nothing, Int64}
        density::Density{<: DensityDependence}
        dependency::Union{Nothing, Type{<:LifeStage}}
        N0::Int64
    end

Data for life stages. Applies to any organism represented by stage-structured population equations.

# Arguments

  - `μ_temperature_response::TemperatureResponse`: Mortality rate. Responsive to temperature.
  - `q_temperature_response::TemperatureResponse`: Developmental rate. 1/total time (days or portion thereof) spent in stage. Responsive to temperature.
  - `n::Union{Nothing, Int64}`: Number of bins allocated to stage (parameter, Erlang distribution).
  - `density::Density`: Specify density dependence model.
  - `dependency::Union{Nothing, Type{<:LifeStage}}`: Organism internal reference, do not modify.
  - `N0::Int64`: Initial stage-specific population count per node. Specify "0" for all stages except `Female` life stage.
"""
mutable struct Stage{L <: LifeStage}
    μ_temperature_response::TemperatureResponse
    q_temperature_response::TemperatureResponse
    n::Union{Nothing, Int64}
    density::Density{<:DensityDependence}
    dependency::Union{Nothing, Type{<:LifeStage}}
    N0::Int64
end

"""
    Stage{L}(μ::Float64, n::Int64, density, dependency, N0::Int64) where {L <: LifeStage}

Return juvenile life stage populated with input data. Applicable to holometabolous (complete) or hemimetabolous (partial) metamorphosing insect species.
"""
function Stage{L}(
    μ::Float64,
    q::Float64,
    n::Int64,
    density,
    dependency,
    N0::Int64,
) where {L <: LifeStage}
    return Stage{L}(NoResponse(μ), NoResponse(q), n, density, dependency, N0)
end

"""
    struct Egg <: LifeStage end

`Egg` life stage. Applicable to holometabolous (complete) or hemimetabolous (partial) metamorphosing insect species.
"""
struct Egg <: LifeStage end

"""
    struct Larva <: LifeStage end

`Larva` life stage. Also referred to as "nymph" life stage. Applicable to holometabolous (complete) or hemimetabolous (partial) metamorphosing insect species.
"""
struct Larva <: LifeStage end

"""
    struct Pupa <: LifeStage end

`Pupa` life stage. Applicable to holometabolous (complete) metamorphosing insect species.
"""
struct Pupa <: LifeStage end

"""
    struct Male <: LifeStage end

`Male` life stage. Applicable to holometabolous (complete) or hemimetabolous (partial) metamorphosing insect species.
"""
struct Male <: LifeStage end

"""
    Stage{Male}(μ::Float64, n::Int64, density, dependency, N0::Int64)

Return `Male` life stage populated with input data. Not dynamically responsive to temperature. Applicable to holometabolous (complete) or hemimetabolous (partial) metamorphosing insect species.
"""
function Stage{Male}(μ::Float64, n::Int64, density, dependency, N0::Int64)
    return Stage{Male}(NoResponse(μ), NoResponse(0.0), n, density, dependency, N0)
end

"""
    Stage{Male}(μ::TemperatureResponse, n::Int64, density, dependency, N0::Int64)

Return `Male` life stage populated with input data. Dynamically responsive to temperature. Applicable to holometabolous (complete) or hemimetabolous (partial) metamorphosing insect species.
"""
function Stage{Male}(μ::TemperatureResponse, n::Int64, density, dependency, N0::Int64)
    return Stage{Male}(μ, NoResponse(0.0), n, density, dependency, N0)
end

"""
    struct Female <: LifeStage end

Female life stage. Applicable to holometabolous (complete) or hemimetabolous (partial) metamorphosing insect species.
"""
struct Female <: LifeStage end

"""
    Stage{Female}(μ::Float64, n::Int64, density, dependency, N0::Int64)

Return `Female` life stage populated with input data. Not dynamically responsive to temperature. Applicable to holometabolous (complete) or hemimetabolous (partial) metamorphosing insect species.
"""
function Stage{Female}(μ::Float64, n::Int64, density, dependency, N0::Int64)
    return Stage{Female}(NoResponse(μ), NoResponse(0.0), n, density, dependency, N0)
end

"""
    Stage{Female}(μ::TemperatureResponse, n::Int64, density, dependency, N0::Int64)

Return `Female` life stage populated with input data. Dynamically responsive to temperature. Applicable to holometabolous (complete) or hemimetabolous (partial) metamorphosing insect species.
"""
function Stage{Female}(μ::TemperatureResponse, n::Int64, density, dependency, N0::Int64)
    return Stage{Female}(μ, NoResponse(0.0), n, density, dependency, N0)
end
