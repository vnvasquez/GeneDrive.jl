module GeneDrive

include("life_stages.jl")
include("genetics.jl")
include("temperature_responses.jl")
include("temperature.jl")

include("data_temperature.jl")

# Structs
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

export NoDensity
export NoResponse
export NoTemperature

export Pupa
export PupaDurationAbiodun
export PupaDurationMoustaid
export PupaDurationRossi

export PupaMortalityAbiodun
export PupaMortalityMoustaid
export PupaMortalityRossi

export SinusoidalTemperature
export Stage

export TimeSeriesTemperature

# Functions
export compute_density

export get_temperature_response
export get_temperature_value



end
