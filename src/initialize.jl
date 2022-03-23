################################################################################
#                                Init Population                               #
################################################################################

"""
        function init_pupa!(node::Node, species::Type{<:Species}, gene_index::Int64, inputs, t)

    Returns initialized pupa life stage.
"""
function init_pupa!(node::Node, species::Type{<:Species},
    gene_index::Int64, inputs, t)

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)

    female = get_lifestage(node, species, Female)
    NF = female.N0
    μF, _ = temperature_effect(ctemp, female)

    pupa = get_lifestage(node, species, Pupa)
    nP = pupa.n
    μP, qP = temperature_effect(ctemp, pupa)

    genetics = get_genetics(node, species)
    gN = count_genotypes(genetics)
    Φ = genetics.Φ

    P0 = zeros(nP, gN)
    P0[end, gene_index] = (NF*μF) / (nP*qP*Φ[gene_index])

    for i in nP-1:-1:1
        P0[i, gene_index] = ((μP + qP*nP)/(qP*nP)) * P0[i+1,gene_index] 
    end

    return P0

end

"""
        function init_egg!(node::Node, species::Type{<:Species}, gene_index::Int64, inputs, t)

    Returns initialized egg life stage.
"""
function init_egg!(node::Node, species::Type{<:Species}, gene_index::Int64, inputs, t)

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)

    female = get_lifestage(node, species, Female)
    NF = female.N0

    egg = get_lifestage(node, species, Egg)
    nE = egg.n
    μE, qE = temperature_effect(ctemp, egg)

    genetics = get_genetics(node, species)
    gN = count_genotypes(genetics)
    Β = genetics.Β

    E0 = zeros(nE, gN)
    E0[1,gene_index] = (Β[gene_index]*NF)/(μE + qE*nE) 

    for i in 2:size(E0)[1]
        E0[i, gene_index] = (qE*nE*E0[i-1, gene_index]) / (μE + qE*nE) 
    end
    return E0

end


"""
        function init_larva!(node::Node, species::Type{<:Species}, gene_index::Int64, inputs, t)

    Returns initialized larva life stage.
"""
function init_larva!(node::Node, species::Type{<:Species}, E0, P0, gene_index::Int64, inputs, t)

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)

    female = get_lifestage(node, species, Female)
    NF = female.N0

    pupa = get_lifestage(node, species, Pupa)
    nP = pupa.n
    μP, qP = temperature_effect(ctemp, pupa)

    egg = get_lifestage(node, species, Egg)
    nE = egg.n
    μE, qE = temperature_effect(ctemp, egg)

    larva = get_lifestage(node, species, Larva)
    nL = larva.n
    μL, qL = temperature_effect(ctemp, larva)

    genetics = get_genetics(node, species)
    gN = count_genotypes(genetics)

    L0 = zeros(nL, gN)
    L0[end, gene_index] = ((μP + qP*nP)/(qL*nL)) * P0[1, gene_index] 

    Lend = L0[end, gene_index]
    Eend = E0[end, gene_index]
    for i in nL-1:-1:1
        L0[i,gene_index] = ((Lend^(i/nL)) * (Eend^((nL-i)/nL)) * (nE^((nL-i)/nL)) * (qE^((nL-i)/nL))) / ((nL^((nL-i)/nL)) * (qL^((nL-i)/nL)))
    end

    return L0

end


"""
        function init_density_dependence!(node::Node, species::Type{<:Species}, eqpop, gene_index::Int64, inputs, t)

    Returns initialized density dependence.
"""
function init_density_dependence!(node::Node, species::Type{<:Species}, eqpop, gene_index::Int64, inputs, t)

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)

    genetics = get_genetics(node, species)
    gN = count_genotypes(genetics)

    egg = get_lifestage(node, species, Egg)
    nE = egg.n
    _, qE = temperature_effect(ctemp, egg)

    larva = get_lifestage(node, species, Larva)
    nL = larva.n
    μL, qL = temperature_effect(ctemp, larva)

    nTot = sum(count_substages(node, species)[1:4]) + count_substages(node, species)[5]*gN
    eqpop = reshape(collect(eqpop),nTot,gN)
    E0 = eqpop[1:nE, 1:gN]
    L0 = eqpop[nE+1:nE+nL, 1:gN]
    Eend = E0[end, gene_index]

    density_param_KL = sum(L0) / ((qE * nE * Eend) / (μL * L0[1,gene_index]) - ((qL * nL) / μL) - 1) 

    density_param_γL = μL / density_param_KL  

    return density_param_KL, density_param_γL

end


"""
        function init_male!(node::Node, species::Type{<:Species}, gene_index::Int64, inputs, t)

    Returns initialized adult male life stage.
"""
function init_male!(node::Node, species::Type{<:Species},
    P0, gene_index::Int64, inputs, t)

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)

    pupa = get_lifestage(node, species, Pupa)
    nP = pupa.n
    μP, qP = temperature_effect(ctemp, pupa)

    male = get_lifestage(node, species, Male)
    μM, _ = temperature_effect(ctemp, male)

    genetics = get_genetics(node, species)
    gN = count_genotypes(genetics)
    Φ = genetics.Φ

    M0 = zeros(1,gN)

    M0[1, gene_index] = ((1-Φ[gene_index])*qP*nP*P0[end, gene_index]) / μM

    return M0

end

"""
        function init_female!(node::Node, species::Type{<:Species}, gene_index::Int64, inputs, t)

    Returns initialized adult female life stage.
"""
function init_female!(node::Node, species::Type{<:Species},
    gene_index::Int64, inputs, t)

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)
    
    female = get_lifestage(node, species, Female)
    F0_begin = female.N0

    genetics = get_genetics(node, species)
    gN = count_genotypes(genetics)

    F0 = zeros(gN, gN)
    F0[gene_index, gene_index] = F0_begin

    return F0

end

"""
        function init_node!(node::Node)

    Returns initialized `Node`.
"""
function init_node!(node::Node)

    t = 0.0
    inputs = ExogenousInputs(node)

    u0 = Vector{Any}(undef, length(keys(node.organisms)))
    densitydep0 = Dict{Any,Any}()

    for (index_organism, key_species) in enumerate(keys(node.organisms))

        genetics = get_genetics(node, key_species)
        gene_index = findfirst(isodd, genetics.all_wildtypes)

        P0 = init_pupa!(node, key_species, gene_index, inputs, t)
        E0 = init_egg!(node, key_species, gene_index, inputs, t)
        L0 = init_larva!(node, key_species, E0, P0, gene_index, inputs, t)
        M0 = init_male!(node, key_species, P0, gene_index, inputs, t)
        F0 = init_female!(node, key_species, gene_index, inputs, t)

        u0_first_guess_perstage = vcat(E0, L0, P0, M0, F0)
        u0[index_organism] = u0_first_guess_perstage

    end

    u0_first_guess_collectedstages = RecursiveArrayTools.ArrayPartition(u0...)

    for (index_organism, key_species) in enumerate(keys(node.organisms))

        genetics = get_genetics(node, key_species)
        gene_index = findfirst(isodd, genetics.all_wildtypes)

        KL, γL = init_density_dependence!(node, key_species, u0_first_guess_collectedstages, gene_index, inputs, t)

        densitydep0[index_organism] = (KL = KL, γL = γL)
        update_density_parameter(node, key_species, Larva; new_param_value = KL)

    end

    eq_pop = NLsolve.nlsolve((du,u) -> population_model_node(du, u, (node, inputs), 0),
    u0_first_guess_collectedstages)
    eqpop = eq_pop.zero

    return eqpop, densitydep0

end

"""
        function init_network!(network::Network)

    Returns initialized `Network`.
"""
function init_network!(network::Network)

    n_nodes = count_nodes(network)
    u0_network = Vector{Any}(undef, n_nodes)
    densitydep0_network = Vector{Any}(undef, n_nodes)

    for (index_node, node_pair) in enumerate(network.nodes)
        key_nodename, node = node_pair
        u0_network[index_node], densitydep0_network[index_node] = init_node!(node)
    end

    u0_first_guess_collectednodes = RecursiveArrayTools.ArrayPartition(u0_network...)
    inputs = ExogenousInputs(network)
    eq_pop = NLsolve.nlsolve((du,u) -> population_model_network(du, u, (network, inputs), 0), u0_first_guess_collectednodes)
    eqpop_network = eq_pop.zero

    return eqpop_network, densitydep0_network

end
