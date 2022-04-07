
```@index
Modules = [GeneDrive]
Pages   = ["features.md"]
```
## [Design](@id features)

`GeneDrive.jl` is predicated on composability. This system design principle allows specific user requirements to be met by mixing and matching modular components (or adding new ones). Because environmental variation - not to mention biological and genetic diversity - can be expressed in many forms, the "composability" of `GeneDrive.jl` means that the features described below may be selected, assembled, and augmented to develop many unique scientific explorations.   

[add figure]

Feature modularity in `GeneDrive.jl` extends to solution methods, more about both jump and diffeq


thanks to the robust selection of solvers supplied by Dig this package builds on the 



## Climate

To enable experimentation with alternative environmental assumptions, three categories or "types" of temperature (°C) inputs are defined in the `GeneDrive.jl` data model. To view them, run:
```julia
julia> ? 
help?> Temperature
```
The `ConstantTemperature` subtype creates a constant environment for the duration of a simulation: `ConstantTemperature(27.0)`. The `SinusoidalTemperature` subtype furnishes an idealized, seasonally fluctuating temperature regime according to user-supplied values for amplitude, periodicity, time period, and mean: `SinusoidalTemperature(4.75, 2, 365, 20.75)`. The `TimeSeriesTemperature` subtype is applicable for empirical data: `TimeSeriesTemperature([vector_of_example_data])`.   

Of course, no treatment of environmental variation would be complete without a feature explicitly acknowledging the realities of climate change. Beyond increases in the average trend that can be accommodated within the three `Temperature` types above, `GeneDrive.jl` accommodates heatwaves and cold snaps via `TemperatureShockData` that imposes time-bound increases or decreases in °C, e.g., augmenting the ambient temperature by 2°C between days 21 and 25 and again between days 37 and 42 in the location of interest: `TemperatureShockData(node,[(21.0, 25.0),(37.0,42.0)],2.0)`.

## Biology

Lifecycle dynamics vary among species thanks to population-specific vital rates and regulatory mechanisms as well as ecological and other factors. The `Organism` type in `GeneDrive.jl` defines the `LifeStages` and [`Genetics`](@ref genetics) particular a user's study species, and is categorized by `Species`. To view the species already defined in the package, run: 
```julia
julia> ? 
help?> Species 
```
`GeneDrive.jl` draws on the thermal biology literature to provide empirically-derived functions that reflect species-appropriate responses to temperature. The following fully parameterized functions - named for the first author of the scientific publication from which they were sourced - are `GeneDrive.jl` data models supplied ready for use: `stages_rossi()`<sup>1</sup> and `stages_moustaid()`<sup>2</sup> are applicable to `AedesAegypti`; `stages_abiodun()`<sup>3</sup> is applicable to `AnophelesGambiae`.

To permit experimentation with alternative field or laboratory-sourced vital rates, the pre-defined `GeneDrive.jl` data model `stages_noresponse()` - named to indicate that no functional response to temperature is assumed - can be populated with stage-specific information using package `accessor` functions. See [API Reference](@ref api). Stage-specific regulatory mechanisms and their contributions to population dynamics can be explored using the various functional forms of density dependence supplied. To view options, run: 
```julia
julia> ? 
help?> DensityDependence 
``` 

## [Genetics](@id genetics)
not just constructs - also talk about supplying new values 

## Geography
network, node, coordinates to overlay with maps/take advantage of cool julia features like LightGraphs 
in the interest of expanding the type of scientific questions that can be asked using basic population models. 

## Interventions


`Release` = both fixed and varied also `ProportionalRelease`


1. insert citation 
2. insert citation 
3. insert citation 




