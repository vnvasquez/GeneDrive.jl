```@index
Modules = [GeneDrive]
Pages   = ["decision_tutorials.md"]
```
# [Decision Model](@id decision_model)

The following shows how to build and solve optimization problems in `GeneDrive.jl`. Decision models allow us to determine the best (optimal) strategy for achieving a goal (objective), taking into account the limitations (constraints) of our system of interest.

To set up the decision model, it is necessary to define an additional component for our data model: our problem's operational limitations (constraints). Operational constraints enter the simulation as fields of the `ReleaseStrategy` struct and are defined per node. View struct fields by running:
```
julia> ? ReleaseStrategy
```
Biological constraints are pre-defined within the decision model. For brevity, we will specify our example operational constraints on top of the `node3` data model created previously.
``` @setup decision_example
using GeneDrive

species = AedesAegypti
genetics = genetics_ridl()
enviro_response = stages_rossi()
update_population_size(enviro_response, 500)
organisms = make_organisms(species, genetics, enviro_response)
temperature = example_temperature_timeseries
coordinates3 = (16.9203, 145.7710)
node3 = Node(:Cairns, organisms, temperature, coordinates3)
release_genotype = get_homozygous_modified(node3, species)
tspan = (1,365);
```

```@example decision_example
# Define constraints using desired fields (re-use `release_genotype`)
example3 = ReleaseStrategy(
    release_this_gene_index = release_genotype,
    release_this_life_stage = Male, 
    release_start_time = 7,
    release_size_max_per_timestep = 1000
); 

# Assign constraints to dict of nodes
my_org_strat = [1 => example3]
my_node_strat = NodeStrategy(1, my_org_strat)

# Build the optimization problem (re-use `node3`, `tspan`)
prob = create_decision_model(
    node3, 
    tspan; 
    node_strategy = my_node_strat
    node_species = [species]
);
```

To solve the decision model as an optimization, a goal (objective) must be supplied. However, even in the absence of an objective function we can derive useful information: without an objective, the solution method auto-selected by `GeneDrive.jl` acts as a nonlinear solver and allows us to compare the behavior of our dynamic and decision models.
```@example decision_example
# Solve
sol = solve_decision_model(prob);

# Format all results for analysis
results = format_decision_model_results(sol);
```
To visualize a subset of the results, run `plot_decision_ridl_females(sol)`. Because no objective function was specified, this output should qualitatively match those from the dynamic model when no intervention is conducted (i.e. when using the same data model, and no RIDL release object is passed to `solve_dynamic_model`).
