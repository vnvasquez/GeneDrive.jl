```@index
Modules = [GeneDrive]
Pages   = ["datasetup_tutorials.md"]
```
# [Data Model](@id data_model)

The first step in any empirical effort is to clean and organize the data. This is also true for computational experiments! `GeneDrive.jl` uses structs to enforce consistency, define relationships, and dynamically assign methods to data. 

Importantly, the struct approach enables modularity: users can construct experiments in a "building block" fashion by assembling information that has already been stored as proper `GeneDrive.jl` inputs. The [Features](@ref features) section outlines the environmental, biological, and genetic details that can be defined thanks to the data model.

Once the information for an experiment has been organized using the data model, we are ready to:
* Save or share our data in a structured and reproducible way.
* Call [Ordinary Differential Equation (ODE) solution methods](@ref dynamic_model) on our data.
* Call [optimization solution methods](@ref decision_model) on our data. 

The code below shows how to construct an example study population using data that is included with the package. 

```@example 
using GeneDrive

# Select species type 
species = AedesAegypti 

# Define how genetic information is passed on 
genetics = genetics_mendelian();

# Choose functional form of environmental response for species life stages
enviro_response = stages_rossi();

# Update population size as desired
update_population_size(enviro_response, 500); 

# Assemble organism
organisms = make_organisms(species, genetics, enviro_response);
```

To fully define an experiment, additional information is relevant: the spatial structure of the population, the ambient temperature of the habitat, and its geographic location should also be defined. The code below demonstrates how to do this; as above, it draws on pre-structured data from `GeneDrive.jl`.

```@example 
# Define temperature functional form and data 
temperature = example_temperature_timeseries;

# Specify the geographic coordinates 
coordinates = (16.1820, 145.7210);

# Define the spatial structure, name the location, and "populate" it 
node1 = Node(:YorkeysKnob, organisms, temperature, coordinates);
```

If the desired spatial structure is a network, we must also define migration rates for subsets of the population that move from node to node within that network. Migration is defined as a nested dictionary wherein the rate at which each genotype and lifestage moves between locations can be optionally specified. When migration rates are not defined for adjacent nodes or specific life stages (e.g., eggs) and genotypes, the default rate is set to zero.
```@example 
# Define a second node 
coordinates2 = (17.0966, 145.7747);
node2 = Node(:Gordonsvale, organisms, temperature, coordinates2);

# Create a network comprised of the two nodes 
network = Network(:Queensland, node1, node2);

# Specify that adult males and females of all genotypes move 
migration_data = Dict( # node1 <-> node2
    ("Male", "AA") => Dict((:YorkeysKnob, :Gordonsvale) => 0.02,
                           (:Gordonsvale, :YorkeysKnob) => 0.02),
    ("Male", "Aa") => Dict((:YorkeysKnob, :Gordonsvale) => 0.02,
                            (:Gordonsvale, :YorkeysKnob) => 0.02),
    ("Male", "aa") => Dict((:YorkeysKnob, :Gordonsvale) => 0.02,
                            (:Gordonsvale, :YorkeysKnob) => 0.02),
    ("Female", "AA") => Dict((:YorkeysKnob, :Gordonsvale) => 0.02,
                            (:Gordonsvale, :YorkeysKnob) => 0.02),
    ("Female", "Aa") => Dict((:YorkeysKnob, :Gordonsvale) => 0.02,
                            (:Gordonsvale, :YorkeysKnob) => 0.02),
    ("Female", "aa") => Dict((:YorkeysKnob, :Gordonsvale) => 0.02,
                            (:Gordonsvale, :YorkeysKnob) => 0.02)
);

# Add migration to the network object 
assign_migration!(network, migration_data, species);
```

For a graphical depiction of the `GeneDrive.jl` data model, see the [Features](@ref features) section of the documentation. 