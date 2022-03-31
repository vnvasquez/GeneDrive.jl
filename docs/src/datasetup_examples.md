# Data Model

```@index
Modules = [GeneDrive]
Pages   = ["datasetup_examples.md"]
```

The first step in any empirical effort is to clean and organize the data. This is also true for computational experiments! `GeneDrive.jl` uses structs to enforce consistency, define relationships, and dynamically assign methods to data. 

Importantly, the struct approach enables modularity: users can construct experiments in a "building block" fashion by assembling information that has already been stored as proper `GeneDrive.jl` inputs. See [XXX data file names XXX] for available data models.   

Once the information for an experiment has been organized using the data model, we are ready to:
* Save or share our data in a structured and reproducible way.
* Call Ordinary Differential Equation (ODE) solution methods on our data. [XXX see dynamic_model.md file XXX]
* Call optimization solution methods on our data. [XXX see decision_model.md file XXX]

The code below shows how to construct an example study population using data that is included with the package. 

```@example 
using DataStructures
using GeneDrive

# Select species type 
species = AedesAegypti 

# Define how genetic information is passed on 
genetics = genetics_mendelian()

# Choose functional form of environmental response for species life stages
enviro_response = stages_rossi()

# Update population size as desired
update_population_size(enviro_response, 500) 

# Assemble organism
organisms = OrderedDict(species => Organism{species}(genetics,enviro_response));
```

To fully define an experiment, additional information is relevant: the spatial structure of the population, the ambient temperature of the habitat, and its geographic location should also be defined. The code below demonstrates how to do this; as above, it draws on pre-structured data from `GeneDrive.jl`.

```@example 
# Define the temperature data as a time series 
temperature = TimeSeriesTemperature(historic1990);

# Specify the geographic coordinates 
coordinates = (16.1820, 145.7210)

# Define the spatial structure, name the location, and "populate" it 
node1 = Node(:YorkeysKnob, organisms, temperature, coordinates)
```

If the desired spatial structure is a network, we must also define migration rates for subsets of the population that move from node to node within that network. Migration is defined as a nested dictionary wherein the rate at which each genotype and lifestage moves between locations can be optionally specified. When migration rates are not defined for adjacent nodes or specific life stages (e.g., eggs) the default is set to zero.
```@example 
# Define a second node 
coordinates2 = (17.0966, 145.7747)
node2 = Node(:Gordonsvale, organisms, temperature, coordinates2)

# Create a network comprised of the two nodes 
network = Network(:Queensland, node1, node2)

# Specify that adult males and females of all genotypes move 
migration_data = Dict(
    # node1 <-> node2
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
)

# Add migration to the network object 
assign_migration!(network2, migration_data, species)
```

For a graphical depiction of the `GeneDrive.jl` data model, see [ XXX figure on Functionalities page XXX ]. 