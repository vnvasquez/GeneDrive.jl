################################################################################
#                               Getters and Setters                            #
################################################################################

########################################
#                 Nodes                #
########################################
function get_nodes(network::Network)
    return network.nodes
end

function count_nodes(network::Network)
    return length(network.nodes)
end

function count_nodes(node::Node)
    return 1
end

########################################
#                Names                 #
########################################

function get_name(node::Node)
    return node.name
end

function get_name(network::Network)
    return network.name
end

########################################
#                Location              #
########################################

function get_location(node::Node)
    return node.location
end

function get_location(network::Network)
    return network.location
end

########################################
#               Organisms              #
########################################

function get_organisms(node::Node)
    return collect(keys(node.organisms))
end

function count_organisms(node::Node)
    return length(node.organisms)
end

function count_organisms(network::Network, node::Symbol)
    return length(network.nodes[node].organisms)
end

function update_organism(node::Node, new_species)
    node.organisms = new_species
    return node
end

function update_organism(network::Network, node::Symbol, new_species)
    network.nodes[node].organisms = new_species
    return network
end

########################################
#              Life Stages             #
########################################

function get_lifestages(node::Node, species::Type{<:Species})
    return node.organisms[species].all_stages
end

function get_lifestage(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})
    return node.organisms[species].all_stages[life_stage]
end

function get_previous_lifestage(node::Node, species::Type{<:Species}, ::Stage{T}) where T <: LifeStage
    prev = node.organisms[species].all_stages[T].dependency
    return node.organisms[species].all_stages[prev]
end

function count_substages(network::Network, node::Symbol, species::Type{<:Species})
    stages_dict = network.nodes[node].organisms[species].all_stages
    # TODO: Improve efficiency
    substage_array = Vector{Int}()
    for (substage, stage) in stages_dict
        substage_array = push!(substage_array,stage.n)
    end
    return substage_array
end

function count_substages(node::Node, species::Type{<:Species})
    stages_dict = node.organisms[species].all_stages
    # TODO: Improve efficiency
    substage_array = Vector{Int}()
    for (substage, stage) in stages_dict
        substage_array = push!(substage_array,stage.n)
    end
    return substage_array
end

function count_substages(stages_dict)
    substage_array = Vector{Int}()
    for (substage, stage) in stages_dict
        substage_array = push!(substage_array,stage.n)
    end
    return substage_array
end

########################################
#               Duration               #
# TODO: update q field to q_temperature_response
########################################


function get_duration(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})
    return node.organisms[species].all_stages[life_stage].q
end

function update_duration(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage}, new_q)
    node.organisms[species].all_stages[life_stage].q = new_q
    return node
end


########################################
#               Mortality              #
# TODO: update μ field to μ_temperature_response
########################################

function get_mortality(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})
    return node.organisms[species].all_stages[life_stage].μ
end

function update_mortality(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage}, new_μ)
    node.organisms[species].all_stages[life_stage].μ = new_μ
    return node
end

########################################
#               Density                #
########################################

function get_density(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})
    return node.organisms[species].all_stages[life_stage].density
end

function update_density_parameter(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage}; new_param_value::Float64)
    dens = node.organisms[species].all_stages[life_stage].density
    dens.param = new_param_value
    return node
end

function update_density_model(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage}; new_density_model::Density)
    dens = node.organisms[species].all_stages[life_stage].density
    dens.model = new_density_model
    return node
end

function update_density(node::Node, species::Type{<:Species}, stage::Type{<:LifeStage}; new_density::Density)
    new_node = deepcopy(node)
    new_node.organisms[species].all_stages[stage].density = new_density
    return new_node
end

########################################
#               Genetics               #
# TODO: update each of these
########################################

function get_genetics(node::Node, species::Type{<:Species})
    return node.organisms[species].gene_data
end

function get_genetics(network::Network, node::Symbol, species::Type{<:Species})
    return network.nodes[node].organisms[species].gene_data
end

function get_genotypes(node::Node, species::Type{<:Species})
    return node.organisms[species].gene_data.all_genotypes
end

function count_genotypes(network::Network, node::Symbol, species::Type{<:Species})
    return length(network.nodes[node].organisms[species].gene_data.all_genotypes)
end

function count_genotypes(node::Node, species::Type{<:Species})
    return length(node.organisms[species].gene_data.all_genotypes)
end

function count_genotypes(genetics::Genetics)
    return length(genetics.all_genotypes)
end

function get_homozygous_modified(node::Node, species::Type{<:Species})
    return findfirst(isodd, node.organisms[species].gene_data.all_modified)
end

function get_homozygous_modified(network::Network, node::Node, species::Type{<:Species})
    return findfirst(isodd, network.nodes[node].organisms[species].gene_data.all_modified)
end

function get_wildtype(node::Node, species::Type{<:Species})
    return findfirst(isodd, node.organisms[species].gene_data.all_wildtypes)
end

function get_wildtype(network::Network, node::Node, species::Type{<:Species})
    return findfirst(isodd, network.nodes[node].organisms[species].gene_data.all_wildtypes)
end

function update_genetics(node::Node, species::Type{<:Species}, new_genetics)
    node.organisms[species].gene_data = new_genetics
    return node
end

function update_genetics_Ω(node::Node, species::Type{<:Species}, new_omega::Array{Float64,1})
    node.organisms[species].gene_data.Ω = new_omega
    return node
end

function update_genetics_Ω(gene_data, new_omega::Array{Float64,1})
    gene_data.Ω = new_omega
    return gene_data
end

function update_genetics_Β(node::Node, species::Type{<:Species}, new_beta::Array{Float64,1})
    node.organisms[species].gene_data.Β = new_beta
    return node
end

function update_genetics_Η(node::Node, species::Type{<:Species}, new_eta::Array{Float64,1})
    node.organisms[species].gene_data.Η = new_eta
    return node
end

function update_genetics_S(node::Node, species::Type{<:Species}, new_sigma::Array{Float64,1})
    node.organisms[species].gene_data.S = new_sigma
    return node
end

########################################
#              Migration               #
########################################

function get_migration(network::Network, node::Symbol, species::Type{<:Species})
    return network.nodes[node].organisms[species].migration
end

function update_migration(network::Network, node::Symbol, species::Type{<:Species}, new_migration)
    network.nodes[node].organisms[species].migration = new_migration
    return network
end

########################################
#             Temperature              #
# TODO: revisit each of these setups
########################################

function get_temperature(node::Node)
    return node.temperature
end

function update_temperature(node::Node, temp_type::Type{<:ConstantTemperature}, new_temperature::Float64)
    node.temperature = ConstantTemperature(new_temperature)
    return node
end

function update_temperature(node::Node, temp_type::Type{<:TimeSeriesTemperature}, new_temperature::Vector{Float64})
    node.temperature = TimeSeriesTemperature(new_temperature)
    return node
end

#= TODO: check the following for update temp:
Fix the sinusoidal function later: requires correct setup to take new_temp as the mean
AND necessary data for the rest of the struct fields

function update_temperature(node::Node, temp_type::Type{<:SinusoidalTemperature}, new_temperature::Float64, Afield, Bfield, Cfield, Dfield)
    new_node = deepcopy(node)
    return new_node.temperature = SinusoidalTemperature(new_temperature)
end
=#
