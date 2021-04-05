module GeneDrive

# Dependencies
#####################
import DataStructures
import DiffEqBase
const diffeq = DiffEqBase
import DiffEqCallbacks
const diffeqCB = DiffEqCallbacks
import NLsolve
import LinearAlgebra
import RecursiveArrayTools

# Files
#####################
# TODO: Random.seed!(123)
include("temperature.jl")
include("temperature_responses.jl")
include("life_stages.jl")
include("genetics.jl")
include("organisms.jl")

#include("temperature.jl")
include("spatial_structure.jl")
include("exogenous_change.jl")

include("accessors.jl")
include("definitions.jl")

include("initialize.jl")
include("population_dynamics.jl")

include("dynamic_model_node.jl")
include("dynamic_model_network.jl")


# Structs
#####################
export AedesAegypti
export AedesAlbopictus
export AnophelesGambiae
export AnophelesArabiensis

export AdultMortalityAbiodun
export AdultMortalityMoustaid
export AdultMortalityRossi

export ConstantTemperature

export Density
export Drive

export Egg
export EggDurationAbiodun
export EggDurationMoustaid
export EggDurationRossi

export EggMortalityAbiodun
export EggMortalityMoustaid
export EggMortalityRossi

export ExogenousInputs

export Female

export Genetics

export Larva
export LarvaDurationAbiodun
export LarvaDurationMoustaid
export LarvaDurationRossi

export LarvaMortalityAbiodun
export LarvaMortalityMoustaid
export LarvaMortalityRossi

export LogisticDensity
export LifeStage
export LinearDensity

export Male

export Network
export Node

export NoDensity
export NoResponse

export Organism

export Pupa
export PupaDurationAbiodun
export PupaDurationMoustaid
export PupaDurationRossi

export PupaMortalityAbiodun
export PupaMortalityMoustaid
export PupaMortalityRossi

export Release

export SinusoidalTemperature
export Species
export Stage

export Temperature
export TemperatureResponse
export TemperatureSeriesData
export TemperatureShockData
export TimeSeriesTemperature

# Genetic structs
#####################
export WW
export ww
export HH
export Hh
export HR
export hh
export hR
export RR

# Functions
#####################
export assign_migration!

export compute_density
export count_genotypes
export count_nodes

export count_organisms
export count_substages

export get_density
export get_duration
export get_exogenous_intervention
export get_exogenous_temperature

export get_genetics
export get_genotypes
export get_homozygous_modified

export get_lifestage
export get_lifestages
export get_location

export get_migration
export get_mortality
export get_name
export get_nodes

export get_organisms

export get_previous_lifestage

export get_temperature
export get_temperature_response
export get_temperature_value

export get_wildtype

export init_node!
export init_network!

export Network

export population_model_node
export population_model_network

export Release

export temperature_effect

export update_density
export update_density_model
export update_density_parameter
export update_duration

export update_genetics
export update_genetics_Β
export update_genetics_Η
export update_genetics_Ω
export update_genetics_S

export update_migration
export update_mortality

export update_organism
export update_temperature

end
