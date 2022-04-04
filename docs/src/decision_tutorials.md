```@index
Modules = [GeneDrive]
Pages   = ["decision_tutorials.md"]
```
# [Decision Model](@id decision_model)

The following shows how to build and solve optimization problems in `GeneDrive.jl`. Decision models allow us to determine the best (optimal) strategy for achieving a goal (objective), taking into account the limitations (constraints) of our system of interest. 

To set up the decision model, it is necessary to define an additional component for our data model: our problem's operational limitations (constraints). Operational constraints enter the simulation as fields of the `ReleaseStrategy` struct and are defined per node. View struct fields by running: 
```julia
julia> ? 
help?> ReleaseStrategy
```

Biological constraints are pre-defined within the decision model. For brevity, we will specify our example operational constraints on top of the `node3` data model created previously. 
```@example
# Define constraints using desired fields (re-use `release_genotype`)
node3_strategy = ReleaseStrategy(release_this_gene_index = release_genotype, 
    release_this_life_stage = Male, release_start_time = 7, 
    release_size_max_per_timestep = 1000)

# Assign constraints to dict of nodes 
mystrategy = Dict(1 => node3_strategy)

# Build the optimization problem (re-use `node3`, `tspan`)
prob = create_decision_model(node3, tspan; node_strategy = mystrategy);
```

To solve the decision model as an optimization, a goal (objective) must be supplied. However, even in the absence of an objective function we can derive useful information: without an objective, the solution method auto-selected by `GeneDrive.jl` acts as a nonlinear solver and allows us to compare the behavior of our dynamic and decision models. 
```@example 
# Solve 
sol = solve_decision_model(prob);

# Format results (use helper function)
results = make_decision_model_traces(sol);

# Visualize
plot(results) 
```

