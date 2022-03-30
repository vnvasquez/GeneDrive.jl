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

import JuMP
import Ipopt
import Gurobi
import Juniper

# Files
#####################
# TODO: Random.seed!(123)
include("temperature.jl")
include("temperature_responses.jl")

include("life_stages.jl")
include("genetics.jl")
include("organisms.jl")

include("spatial_structure.jl")
include("exogenous_change.jl")

include("accessors.jl")
include("definitions.jl")

include("initialize.jl")
include("population_dynamics.jl")

include("dynamic_model_node.jl")
include("dynamic_model_network.jl")

include("constraints.jl")
include("decision_model.jl")

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
export ProportionalRelease
export ReleaseStrategy

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
export Wolbachia
export WW
export ww

export RIDL
export MCR
export HH
export Hh
export HR
export hh
export hR
export RR

export SLHoming
export WH
export WR
export WB
export HB
export RB
export BB

export SplitDrive
export WWWW
export WWWH
export WWWR
export WWWB
export WWHH
export WWHR
export WWHB
export WWRR
export WWRB
export WWBB
export WCWW
export WCWH
export WCWR
export WCWB
export WCHH
export WCHR
export WCHB
export WCRR
export WCRB
export WCBB
export CCWW
export CCWH
export CCWR
export CCWB
export CCHH
export CCHR
export CCHB
export CCRR
export CCRB
export CCBB

# Functions
#####################
export assign_migration!

export compute_density
export count_genotypes
export count_nodes

export count_organisms
export count_substages

export create_decision_model

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
export get_names
export get_nodes

export get_organisms

export get_previous_lifestage

export get_initial_temperature
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

export solve_decision_model

export temperature_effect

export update_density!
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
export update_population_size
export update_temperature

end
