################################################################################
#                               Getters and Setters                            #
################################################################################

########################################
#                 Nodes                #
########################################
"""
    get_nodes(network::Network)

Return the `Node` object(s) contained in the `Network` object.
"""
function get_nodes(network::Network)
    return network.nodes
end

"""
    count_nodes(network::Network)

Return a count of the nodes contained in the `Network` object.
"""
function count_nodes(network::Network)
    return length(network.nodes)
end

"""
    count_nodes(node::Node)

Return 1 (the count of the nodes contained in the `Node` object).
"""
function count_nodes(node::Node)
    return 1
end

########################################
#                Names                 #
########################################

"""
    get_name(node::Node)

Return the name (symbol) of the `Node` object.
"""
function get_name(node::Node)
    return node.name
end

"""
    get_name(network::Network)

Return the name (symbol) of the `Network` object.
"""
function get_name(network::Network)
    return network.name
end

"""
    get_names(network::Network)

Return the names (symbols) of the nodes within the `Network` object.
"""
function get_names(network::Network)
    return keys(network.nodes)
end

########################################
#                Location              #
########################################

"""
    get_location(node::Node)

Return the coordinates of `Node`.
"""
function get_location(node::Node)
    return node.location
end

"""
    get_location(network::Network)

Return the coordinates of `Network`.
"""
function get_location(network::Network)
    return network.location
end

########################################
#               Organisms              #
########################################

"""
    get_organisms(node::Node)

Return vector of species populating `Node`.
"""
function get_organisms(node::Node)
    return collect(keys(node.organisms))
end
"""
    count_organisms(node::Node)

Count the species populating `Node`.
"""
function count_organisms(node::Node)
    return length(node.organisms)
end

"""
    get_organisms(network::Network, node::Symbol)

Return vector of species populating the specified `Node` of `Network`.
"""
function get_organisms(network::Network, node::Symbol)
    return collect(keys(network.nodes[node].organisms))
end

"""
    count_organisms(network::Network, node::Symbol)

Count the species populating the specified `Node` of `Network`.
"""
function count_organisms(network::Network, node::Symbol)
    return length(network.nodes[node].organisms)
end

"""
    update_organism(node::Node, new_species)

Update the species populating `Node`.
"""
function update_organism(node::Node, new_species)
    node.organisms = new_species
    return node
end

"""
    update_organism(network::Network, node::Symbol, new_species)

Update the species populating the specified `Node` of `Network`.
"""
function update_organism(network::Network, node::Symbol, new_species)
    network.nodes[node].organisms = new_species
    return network
end

########################################
#              Life Stages             #
########################################
"""
    get_lifestages(node::Node, species::Type{<:Species})

Return all data for all `LifeStage`s of `Species` in `Node`.
"""
function get_lifestages(node::Node, species::Type{<:Species})
    return node.organisms[species].all_stages
end

"""
    get_lifestage(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})

Return data for specified `LifeStage` of `Species` in `Node`.
"""
function get_lifestage(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})
    return node.organisms[species].all_stages[life_stage]
end

"""
    get_previous_lifestage(node::Node, species::Type{<:Species}, ::Stage{T}) where T <: LifeStage

Show `LifeStage` dependency: return data for the `LifeStage` previous to the specified `LifeStage` of `Species` in `Node`.
"""
function get_previous_lifestage(
    node::Node,
    species::Type{<:Species},
    ::Stage{T},
) where {T <: LifeStage}
    prev = node.organisms[species].all_stages[T].dependency
    return node.organisms[species].all_stages[prev]
end

"""
    count_substages(network::Network, node::Symbol, species::Type{<:Species})

Return count of the total substages for `Species` in the specified `Node` of `Network`.
"""
function count_substages(network::Network, node::Symbol, species::Type{<:Species})
    stages_dict = network.nodes[node].organisms[species].all_stages
    # TODO: Improve efficiency
    substage_array = Vector{Int}()
    for (substage, stage) in stages_dict
        substage_array = push!(substage_array, stage.n)
    end
    return substage_array
end

"""
    count_substages(node::Node, species::Type{<:Species})

Return count of the total substages for `Species` in `Node`.
"""
function count_substages(node::Node, species::Type{<:Species})
    stages_dict = node.organisms[species].all_stages
    # TODO: Improve efficiency
    substage_array = Vector{Int}()
    for (substage, stage) in stages_dict
        substage_array = push!(substage_array, stage.n)
    end
    return substage_array
end

"""
    count_substages(node::Node, species::Type{<:Species}, stage::Type{<:LifeStage})

Return count of the substages in `LifeStage` for `Species` in `Node`.
"""
function count_substages(node::Node, species::Type{<:Species}, stage::Type{<:LifeStage})
    substages = node.organisms[species].all_stages[stage].n
end

"""
    count_substages(node::Node, species::Type{<:Species}, stage::Type{Female})

Return count of the substages for `Species` in `Node`. Specific to the `Female` lifestage.
"""
function count_substages(node::Node, species::Type{<:Species}, stage::Type{Female})
    substages = node.organisms[species].all_stages[stage].n
    gN = count_genotypes(node, species)
    substages_total = substages * gN
end

function count_substages(stages_dict)
    substage_array = Vector{Int}()
    for (substage, stage) in stages_dict
        substage_array = push!(substage_array, stage.n)
    end
    return substage_array
end

"""
    update_population_size(stages::OrderedDict{DataType, Stage}, new_popsize::Int64)

Return updated population size. Note: `new_popsize` argument refers specifically to the Female population; if e.g. new_popsize = 500, the full adult population (Females and Males) will be 500*2.
"""
function update_population_size(stages::DataStructures.OrderedDict{DataType, Stage}, new_popsize::Int64)
    stages[Female].N0 = new_popsize
    return stages
end

########################################
#               Duration               #
########################################

"""
    get_duration(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})

Return the temperature-sensitive duration (`q_temperature_response`) for the `LifeStage` of `Species` in `Node.`
"""
function get_duration(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})
    return node.organisms[species].all_stages[life_stage].q_temperature_response
end

"""
    update_duration(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage}, new_q)

Update the temperature-sensitive duration (`q_temperature_response`) for the `LifeStage` of `Species` in `Node.`
"""
function update_duration(
    node::Node,
    species::Type{<:Species},
    life_stage::Type{<:LifeStage},
    new_q,
)
    node.organisms[species].all_stages[life_stage].q_temperature_response = new_q
    return node
end

########################################
#               Mortality              #
########################################

"""
    get_mortality(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})

Return the temperature-sensitive mortality (`μ_temperature_response`) for the `LifeStage` of `Species` in `Node.`
"""
function get_mortality(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})
    return node.organisms[species].all_stages[life_stage].μ_temperature_response
end

"""
    update_mortality(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage}, new_μ)

Update the temperature-sensitive mortality (`μ_temperature_response`) for the `LifeStage` of `Species` in `Node.`
"""
function update_mortality(
    node::Node,
    species::Type{<:Species},
    life_stage::Type{<:LifeStage},
    new_μ,
)
    node.organisms[species].all_stages[life_stage].μ_temperature_response = new_μ
    return node
end

########################################
#               Density                #
########################################

"""
    get_density(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})

Return the density dependence model and parameterization for `LifeStage` of `Species` in `Node`.
"""
function get_density(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})
    return node.organisms[species].all_stages[life_stage].density
end

"""
    update_density_parameter(node::Node, species::Type{<:Species}, ::Type{T}; new_param_value::Float64) where T <: LifeStage

Update the parameterization of the density dependence model for `LifeStage` of `Species` in `Node`.
"""
function update_density_parameter(
    node::Node,
    species::Type{<:Species},
    ::Type{T};
    new_param_value::Float64,
) where {T <: LifeStage}
    @debug "changed the density parameters $T, $new_param_value"
    node.organisms[species].all_stages[T].density.param = new_param_value
    return
end

"""
    update_density_model(node::Node, species::Type{<:Species}, ::Type{T}; new_density_model::Density) where T <: LifeStage

Update the functional form of the density dependence model for `LifeStage` of `Species` in `Node`.
"""
function update_density_model(
    node::Node,
    species::Type{<:Species},
    ::Type{T};
    new_density_model::Density,
) where {T <: LifeStage}
    node.organisms[species].all_stages[T].density.model = new_density_model
    return
end

"""
    update_density!(node::Node, species::Type{<:Species}, ::Type{T}; new_density::Density) where T <: LifeStage

Update both the functional form and parameterization of the density dependence model for `LifeStage` of `Species` in `Node`.
"""
function update_density!(
    node::Node,
    species::Type{<:Species},
    ::Type{T};
    new_density::Density,
) where {T <: LifeStage}
    new_node = deepcopy(node)
    new_node.organisms[species].all_stages[T].density = new_density
    return new_node
end

########################################
#               Genetics               #
# TODO: update each of these
########################################

"""
    get_genetics(node::Node, species::Type{<:Species})

Return `Genetics` data for `Species` in `Node`.
"""
function get_genetics(node::Node, species::Type{<:Species})
    return node.organisms[species].gene_data
end

"""
    get_genetics(network::Network, node::Symbol, species::Type{<:Species})

Return `Genetics` data for `Species` in `Node` of `Network`.
"""
function get_genetics(network::Network, node::Symbol, species::Type{<:Species})
    return network.nodes[node].organisms[species].gene_data
end

"""
    update_genetics(node::Node, species::Type{<:Species}, new_genetics)

Update `Genetics` data for `Species` in `Node`.
"""
function update_genetics(node::Node, species::Type{<:Species}, new_genetics)
    node.organisms[species].gene_data = new_genetics
    return node
end

"""
    update_genetics_Ω(node::Node, species::Type{<:Species}, new_omega::Array{Float64,1})

Update `Ω` data in `Genetics` for `Species` in `Node`.
"""
function update_genetics_Ω(
    node::Node,
    species::Type{<:Species},
    new_omega::Array{Float64, 1},
)
    node.organisms[species].gene_data.Ω = new_omega
    return node
end

"""
    update_genetics_Ω(gene_data, new_omega::Array{Float64,1})

Update `Ω` data in `Genetics`. Helper function.
"""
function update_genetics_Ω(gene_data, new_omega::Array{Float64, 1})
    gene_data.Ω = new_omega
    return gene_data
end

"""
    update_genetics_Β(node::Node, species::Type{<:Species}, new_beta::Array{Float64,1})

Update `Β` data in `Genetics` for `Species` in `Node`.
"""
function update_genetics_Β(
    node::Node,
    species::Type{<:Species},
    new_beta::Array{Float64, 1},
)
    node.organisms[species].gene_data.Β = new_beta
    return node
end

"""
    update_genetics_Η(node::Node, species::Type{<:Species}, new_eta::Array{Float64,1})

Update `Η` data in `Genetics` for `Species` in `Node`.
"""
function update_genetics_Η(node::Node, species::Type{<:Species}, new_eta::Array{Float64, 1})
    node.organisms[species].gene_data.Η = new_eta
    return node
end

"""
    update_genetics_Η(gene_data, new_eta::Array{Float64,1})

Update `Η` data in `Genetics`. Helper function.
"""
function update_genetics_Η(gene_data, new_eta::Array{Float64, 1})
    gene_data.Η = new_eta
    return gene_data
end

"""
    update_genetics_S(node::Node, species::Type{<:Species}, new_sigma::Array{Float64,1})

Update `S` data in `Genetics` for `Species` in `Node`.
"""
function update_genetics_S(
    node::Node,
    species::Type{<:Species},
    new_sigma::Array{Float64, 1},
)
    node.organisms[species].gene_data.S = new_sigma
    return node
end

"""
    get_genotypes(node::Node, species::Type{<:Species})

Return genotype data from `Genetics` for `Species` in `Node`.
"""
function get_genotypes(node::Node, species::Type{<:Species})
    return node.organisms[species].gene_data.all_genotypes
end

"""
    count_genotypes(node::Node, species::Type{<:Species})

Return the total count of genotypes from `Genetics` for `Species` in `Node`.
"""
function count_genotypes(node::Node, species::Type{<:Species})
    return length(node.organisms[species].gene_data.all_genotypes)
end

"""
    count_genotypes(network::Network, node::Symbol, species::Type{<:Species})

Return the total count of genotypes from `Genetics` for `Species` in `Node` of `Network`.
"""
function count_genotypes(network::Network, node::Symbol, species::Type{<:Species})
    return length(network.nodes[node].organisms[species].gene_data.all_genotypes)
end

"""
    count_genotypes(genetics::Genetics)

Return the total count of genotypes in the `Genetics` object.
"""
function count_genotypes(genetics::Genetics)
    return length(genetics.all_genotypes)
end

"""
    get_homozygous_modified(node::Node, species::Type{<:Species})

Return the index of the homozygous modified genotype in `Genetics` for `Species` in `Node`.
"""
function get_homozygous_modified(node::Node, species::Type{<:Species})
    return findfirst(isodd, node.organisms[species].gene_data.all_modified)
end

"""
    get_homozygous_modified(network::Network, node::Node, species::Type{<:Species})

Return the index of the homozygous modified genotype in `Genetics` for `Species` in `Node` of `Network`.
"""
function get_homozygous_modified(network::Network, node::Node, species::Type{<:Species})
    return findfirst(isodd, network.nodes[node].organisms[species].gene_data.all_modified)
end

"""
    get_wildtype(node::Node, species::Type{<:Species})

Return the index of the wildtype genotype in `Genetics` for `Species` in `Node`.
"""
function get_wildtype(node::Node, species::Type{<:Species})
    return findfirst(isodd, node.organisms[species].gene_data.all_wildtypes)
end

"""
    get_wildtype(network::Network, node::Node, species::Type{<:Species})

Return the index of the wildtype genotype in `Genetics` for `Species` in `Node` of `Network`.
"""
function get_wildtype(network::Network, node::Node, species::Type{<:Species})
    return findfirst(isodd, network.nodes[node].organisms[species].gene_data.all_wildtypes)
end

########################################
#              Migration               #
########################################

"""
    get_migration(network::Network, species::Type{<:Species})

Return the migration characterizing each genotype and `Lifestage` for `Species` in `Network`.
"""
function get_migration(network::Network, species::Type{<:Species})
    return network.migration[species]
end

"""
    update_migration(network::Network, species::Type{<:Species}, new_migration)

Update the migration characterizing each genotype and `Lifestage` for `Species` in `Network`.
"""
function update_migration(network::Network, species::Type{<:Species}, new_migration)
    network.migration[species] = new_migration
    return network
end

########################################
#             Temperature              #
########################################

"""
    get_temperature(node::Node)

Return `Temperature` for `Node`. Object contains `type` and `values`.
"""
function get_temperature(node::Node)
    return node.temperature
end

"""
    get_initial_temperature(node::Node)

Return entry for the first index of `values` in `Temperature` of `Node.`
"""
function get_initial_temperature(node::Node)
    ctemp = node.temperature
    initial_ctemp = ctemp.values[1]
    return initial_ctemp
end

"""
    update_temperature(node::Node, temp_type::Type{<:ConstantTemperature}, new_temperature::Float64)

Update the `type` and `values` of `Temperature` for `Node`.
"""
function update_temperature(
    node::Node,
    temp_type::Type{<:ConstantTemperature},
    new_temperature::Float64,
)
    node.temperature = ConstantTemperature(new_temperature)
    return node
end

"""
    update_temperature(node::Node, temp_type::Type{<:TimeSeriesTemperature}, new_temperature::Vector{Float64})

Update the `type` and `values` of `Temperature` for `Node`.
"""
function update_temperature(
    node::Node,
    temp_type::Type{<:TimeSeriesTemperature},
    new_temperature::Vector{Float64},
)
    node.temperature = TimeSeriesTemperature(new_temperature)
    return node
end
