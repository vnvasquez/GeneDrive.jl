################################################################################
#                                     Node                                     #
################################################################################

"""
        mutable struct Node
            name::Symbol
            organisms::DataStructures.OrderedDict{Type{<:Species}, Organism}
            temperature::Temperature
            location::Tuple{Float64,Float64}
        end

    Data for a single spatial node.

# Arguments
- `name::Symbol`: Name of node, usually location-relevant.
- `organisms::DataStructures.OrderedDict{Type{<:Species}, Organism} `: Dictionary containing data for all species inhabiting node.
- `temperature::Temperature`: Climatic specification for temperature.
- `location::Tuple{Float64,Float64}`: Geographic location denoted by coordinates.
"""
mutable struct Node
    name::Symbol
    organisms::DataStructures.OrderedDict{Type{<:Species}, Organism}
    temperature::Temperature
    location::Tuple{Float64,Float64}
end

################################################################################
#                           Network and Migration                              #
################################################################################
"""
        struct Network
            name::Symbol
            nodes::DataStructures.OrderedDict{Symbol,Node}
            migration::DataStructures.OrderedDict{DataType},Array{Matrix{Float64},2}}
            locations_key_map
        end

    Data for multiple interconnected spatial nodes.

# Arguments
- `name::Symbol`: Name of network, usually location-relevant.
- `nodes::DataStructures.OrderedDict{Symbol,Node}`: Dictionary of data for nodes included in network (metapopulation).
- `migration::DataStructures.OrderedDict{DataType}, Array{Matrix{Float64},2}}`: Data defining species, genotype, and stage-specific transition rates. Entries default to zero.
- `locations_key_map`: Mapping of node location information. Used to assign transition rates.
"""
struct Network
    name::Symbol
    nodes::DataStructures.OrderedDict{Symbol,Node}
    migration::DataStructures.OrderedDict{DataType, Array{Matrix{Float64},2}}
    locations_key_map
end

function _make_migration_array(network_dict::AbstractDict)
    migration = DataStructures.OrderedDict{DataType, Array{Matrix{Float64},2}}()
    node_count = length(keys(network_dict))
    for (key_nodename, index_node) in network_dict
        for (key_species, index_organisms) in index_node.organisms
            gene_count = count_genotypes(index_organisms.gene_data)
            stages = index_organisms.all_stages
            stage_count = sum(count_substages(stages)[1:4]) + count_substages(stages)[5] * gene_count
            trans_mat = Matrix{Matrix{Float64}}(undef, stage_count, gene_count)
            for i in eachindex(trans_mat)
                trans_mat[i] = zeros(node_count, node_count)
            end
            migration[key_species] = trans_mat
        end
    end
    return migration
end

"""
        Network(name::Symbol, nodes::DataStructures.OrderedDict{Symbol,Node}...)

    Returns `Network` instance containing specified metapopulation information.
"""
function Network(name::Symbol, nodes::Node...)
    network_dict = DataStructures.OrderedDict{Symbol,Node}()
    for node in nodes
        key_nodename = get_name(node)
        network_dict[key_nodename] = node
    end
    migration = _make_migration_array(network_dict)
    locations_key_map = Dict(key => ix for (ix, key) in enumerate(keys(network_dict)))
    return Network(name, network_dict, migration, locations_key_map)
end

function _fill_migration_array!(initial_migration, migration_data,
    gene_to_index_migration_matrix,
    stage_to_index_migration_matrix,
    locations_key_map)

    for (life_gene_key, migration_matrix) in migration_data
    stage, gene = life_gene_key
    stage_index = stage_to_index_migration_matrix[life_stage_key_map[stage]]
    gene_index = gene_to_index_migration_matrix[genetics_key_map[gene]]

        for ix in stage_index
            mat = initial_migration[ix, gene_index]

            for (from_to_nodes, move_rate) in migration_matrix
                from_node, to_node = from_to_nodes
                from_index = locations_key_map[from_node]
                to_index = locations_key_map[to_node]
                mat[from_index, to_index] = move_rate
            end

            for i in 1:size(mat)[1]
                mat[i, i] = -sum(mat[i, :])
            end
        end
    end
end

"""
        assign_migration!(network::Network, migration_data::Dict, species::Type{<:Species})

- Returns `Network` instance with `migration` field populated by user-specified transition rates.
- Input data must be formatted as a nested dictionary. First level denotes relevant life stage and gene, second level includes to/from nodes and transition rate.
- Stage and gene combinations not specified by input data retain default transition rate of zero.
- See `data_migration.jl` file for example.
"""
function assign_migration!(network::Network, migration_data::Dict, species::Type{<:Species})

    node = first(values(get_nodes(network)))

    genetics = get_genetics(node, species)
    gN = count_genotypes(genetics)
    gene_to_index_migration_matrix = Dict(g.genotype => ix for (ix, g) in enumerate(genetics.all_genotypes))

    n = count_substages(node, species)
    nE = n[1]
    nL = n[2]
    nP = n[3]
    nJuv = nE+nL+nP
    nM = n[4]
    nF = n[5]*gN + nM + nJuv
    stage_to_index_migration_matrix = Dict(Egg => 1:nE,
                                          Larva => 1+nE:nE+nL,
                                          Pupa => 1+nE+nL:nE+nL+nP,
                                          Male => 1+nJuv,
                                          Female => 1+nJuv+nM:nF)

    initial_migration = network.migration[species]

    _fill_migration_array!(initial_migration, migration_data,
                           gene_to_index_migration_matrix, stage_to_index_migration_matrix,
                           network.locations_key_map)
    return
end
