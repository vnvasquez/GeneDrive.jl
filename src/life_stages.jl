
################################################################################
#                              Density Dependence                              #
################################################################################

abstract type DensityDependence end

"""
        mutable struct Density{D <: DensityDependence}
            model::Type{D}
            param::Float64
        end

    Data for density dependence. Applies to all species, all stages.

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

    Density dependence model of type `Logistic Density`.
"""
struct LogisticDensity <: DensityDependence end

"""
        compute_density(data::Density{LogisticDensity}, stage)

    Returns the effect of logistic density dependence on life stage instance.
    Calculated each time step.
"""
function compute_density(data::Density{LogisticDensity}, stage)
    return (1 + (sum(stage)/data.param))
    # TODO: note to Murat that nonlinearity not an issue bc sum is divided by a single parameter
end

"""
        struct LinearDensity <: DensityDependence end

    Density dependence model of type `Linear Density`.
"""
struct LinearDensity <: DensityDependence end

"""
        compute_density(data::Density{LinearDensity}, stage)

    Returns the effect of linear density dependence on life stage instance.
    Calculated each time step.
"""
function compute_density(data::Density{LinearDensity}, stage)
    return data.param * sum(stage)
    # TODO: check against Hancock paper.
end

"""
        struct NoDensity <: DensityDependence end

    Density dependence model of type `No Density`.
"""
struct NoDensity <: DensityDependence end

"""
        compute_density(data::Density{NoDensity}, stage)

    Returns the effect of no density dependence on life stage instance.
    Calculated each time step.
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
        q::Union{Nothing, Float64}
        n::Union{Nothing, Int64}
        μ::Float64
        density::Density
        N0::Int64
        dependency::Union{Nothing, Type{<:LifeStage}}
    end

Data for life stages. Applies to any organism represented by stage-structured population equations.

# Arguments
- `q::Union{Nothing, Float64}`: Rate of change (1/duration). Duration: total time (days) in each substage.
- `n::Union{Nothing, Int64}`: Number of substages (shape parameter, Erlang distribution).
- `μ::Float64`: Mortality rate
- `density::Density`: Specify density dependence model. Default: `LogisticDensity` in `Larva` life stage, `NoDensity` in all other life stages.
- `N0::Int64`: Initial stage-specific population count per node. Specify "0" for all stages except `Female` life stage.
- `dependency::LifeStage`: Organism internal reference, do not modify.
"""
mutable struct Stage{L <: LifeStage}
    q::Union{Nothing, Float64}
    n::Union{Nothing, Int64} #TODO: Define Erlang aspect more thoroughly in docstring.
    μ::Float64
    density::Density
    N0::Int64
    dependency::Union{Nothing, Type{<:LifeStage}}
end

"""
        struct Egg <: LifeStage end

    Egg life stage. Applicable to holometabolous (complete) or hemimetabolous (partial) metamorphosing insect species.
"""
struct Egg <: LifeStage end

"""
        struct Larva <: LifeStage end

    Larva life stage. Also referred to as "nymph" life stage. Applicable to holometabolous (complete) or hemimetabolous (partial) metamorphosing insect species.
"""
struct Larva <: LifeStage end

"""
        struct Pupae <: LifeStage end

    Pupa life stage. Applicable to holometabolous (complete) metamorphosing insect species.
"""
struct Pupa <: LifeStage end

"""
        struct Male <: LifeStage end

    Male life stage. Applicable to holometabolous (complete) or hemimetabolous (partial) metamorphosing insect species.
"""
struct Male <: LifeStage end

"""
        struct Female <: LifeStage end

    Female life stage. Applicable to holometabolous (complete) or hemimetabolous (partial) metamorphosing insect species.
"""
struct Female <: LifeStage end
