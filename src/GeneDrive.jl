module GeneDrive

# Dependencies
import DataStructures

# Files
include("life_stages.jl")
include("genetics.jl")
include("organisms.jl")

#include("definitions.jl")
include("temperature.jl")
include("spatial_structure.jl")

include("temperature_responses.jl")

include("data_temperature.jl")
include("data_migration.jl")

# Structs
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
export LinearDensity

export Male

export Network
export Node

export NoDensity
export NoResponse
export NoTemperature

export Organism

export Pupa
export PupaDurationAbiodun
export PupaDurationMoustaid
export PupaDurationRossi

export PupaMortalityAbiodun
export PupaMortalityMoustaid
export PupaMortalityRossi

export SinusoidalTemperature
export Stage

export Temperature
export TimeSeriesTemperature

# Functions
export assign_migration!

export compute_density

export get_temperature_response
export get_temperature_value

export Network

end
