
```@index
Modules = [GeneDrive]
Pages   = ["features.md"]
```
## [Design](@id features)

`GeneDrive.jl` is predicated on composability. This system design principle allows specific user requirements to be met by mixing and matching modular components (or adding new ones). Because environmental variation - not to mention biological and genetic diversity - can be expressed in many forms, the "composability" of `GeneDrive.jl` means that the features described below may be selected, assembled, and augmented to develop many unique scientific explorations.   

Feature modularity in `GeneDrive.jl` extends to solution methods: the dynamic model is built on top of [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/), enabling users to experiment with the robust suite of solvers in that package. The decision model employs [`JuMP.jl`](https://jump.dev/JuMP.jl/stable/) and can be run using a plethora of free as well as paid solution algorithms for nonlinear and mixed integer problems. Once built, the same data model can be evaluated using either ODE or optimization approaches with no further customization required.  

## Climate

To enable experimentation with alternative environmental assumptions, three categories or "subtypes" of `Temperature` (°C) inputs are defined in the `GeneDrive.jl` data model: 
* `ConstantTemperature` creates a constant environment for the duration of a simulation: `ConstantTemperature(27.0)`. 
* `SinusoidalTemperature` furnishes an idealized, seasonally fluctuating temperature regime according to user-supplied values for amplitude, periodicity, time period, and mean: `SinusoidalTemperature(4.75, 2, 365, 20.75)`. 
* `TimeSeriesTemperature` is applicable for empirical data: `TimeSeriesTemperature([vector_of_example_data])`.   

Beyond increases in the average trend that can be accommodated within the three `Temperature` subtypes above, `GeneDrive.jl` accommodates heatwaves and cold snaps via `TemperatureShockData` that imposes time-bound increases or decreases in °C, e.g., augmenting the ambient temperature by 2°C between days 21 and 25 and again between days 37 and 42 in the location of interest: `TemperatureShockData(node,[(21.0, 25.0),(37.0,42.0)],2.0)`.

## [Biology](@id biology)

Lifecycle dynamics vary among species thanks to population-specific vital rates and regulatory mechanisms as well as ecological and other factors. The `Organism` type in `GeneDrive.jl` is defined by the `LifeStages` and [`Genetics`](@ref genetics) particular to a given study species, and is categorized by `Species`. To view the species already defined in the package, run: 
```julia
julia> ? 
help?> Species 
```
`GeneDrive.jl` draws on the thermal biology literature to provide empirically-derived functions that characterize species-appropriate temperature responses. The following - named for the first author of the scientific publication from which they were sourced - are fully parameterized `GeneDrive.jl` data models: 
* `stages_rossi()` and `stages_moustaid()`, applicable to `AedesAegypti`
* `stages_abiodun()`, applicable to `AnophelesGambiae`

To permit experimentation with alternative field or laboratory-sourced vital rates, the pre-defined `GeneDrive.jl` data model `stages_noresponse()` - named to indicate that no functional response to temperature is assumed - can be populated with stage-specific information using package `accessor` functions. See [API Reference](@ref api). 

Stage-specific regulatory mechanisms and their contributions to population dynamics can also be explored. See `Density` models under the type `DensityDependence`.

## [Genetics](@id genetics)

The `Genetics` type in `GeneDrive.jl` defines the likelihood with which an organism will produce offspring and pass on its genetic material. This likelihood is currently assumed to be equivalent for all wildtype organisms, following [Mendelian](https://en.wikipedia.org/wiki/Mendelian_inheritance) inheritance principles. However, [genetic-based methods](https://fnih.org/what-we-do/geneconvene/about/genetic-biocontrol) of biocontrol - a category that ranges from actual [modification](https://www.nature.com/articles/d41586-019-02087-5) to bacterial [infection](http://www.eliminatedengue.com/our-research/Wolbachia) - alter natural inheritance patterns, with affected genotypes having a greater or lesser ability to propagate. 

Beyond this, `Genetics` data includes information unique to genetic constructs being developed for biocontrol including fertility rates, sex ratios, and varied fitness costs (e.g. how mortality rates for modified organisms are biased with respect to their wild counterparts). Together with the vital rates defined in `LifeStages`, `Genetics` data is used to create the `Organism`s that populate a `GeneDrive.jl` experiment. The following are fully parameterized data models included in the package: 
* `genetics_mendelian()`: For wildtypes. 
* `genetics_ridl()`: For the biocontrol technology RIDL (suppression), developed by Oxitec. 
* `genetics_wolbachia()`: For two biocontrol approaches (i.e., suppression and replacement) that infect organisms with *Wolbachia* (multiple existing strains).
* `genetics_mcr()`: For a gene drive biocontrol technology (replacement) tested exclusively in laboratory environnments. 


## Geography

Connectivity between geographic locations informs population dynamics by allowing organisms to migrate and spatially propagate genetic material. Study species defined by the `GeneDrive.jl` data model may inhabit a `Node` - a single homogenous habitat - or a heterogeneous `Network`. The `Network` is a collection of interconnected `Node`s, each of which are geolocated by coordinates and may exhibit environmental differences. See the data model [tutorial](@ref data_model) for an example treatment of the fields below: 

```@setup geography
using GeneDrive
species = AedesAegypti 
genetics = genetics_mendelian();
enviro_response = stages_rossi();
update_population_size(enviro_response, 500); 
organisms = make_organisms(species, genetics, enviro_response);
temperature = example_temperature_timeseries; 
coordinates = (16.1820, 145.7210);
node1 = Node(:YorkeysKnob, organisms, temperature, coordinates);
```
```@example geography
node1
```
Movement between all `Node` locations in a `Network` is uniquely specified for each species, life stage, and genotype using [`assign_migration!`](@ref migration). This enables users to account for a diversity of demographic and genetic migration tendencies as well as exogenous factors, facilitating scientific investigation in an [important area](https://www.science.org/content/article/windborne-mosquitoes-may-carry-malaria-hundreds-kilometers). 

## Interventions




