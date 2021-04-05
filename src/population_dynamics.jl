################################################################################
#                         Population: ODE Formulation                          #
################################################################################

"""
        function oviposit(F, cube, Τ, S, Β, gene_index::Int64)

    Returns oviposited eggs (count).
"""
function oviposit(F, cube, Τ, S, Β, gene_index::Int64)

    ΒS = Matrix{Float64}(undef, length(Β), length(S))
    for i in 1:length(Β)
        ΒS[i,:] = Β[i].*S'
    end
    # TODO: add genetic parameters necessary for additional constructs
    O = cube[:,:,gene_index].*Τ[:,:,gene_index].*ΒS.*F

    return sum(O)
end


"""
        function create_egg!(dE, E, node::Node, species::Type{<:Species}, eggsnew, gene_index::Int64, inputs::ExogenousInputs, t)

    Returns egg state.
"""
function create_egg!(dE, E, node::Node, species::Type{<:Species}, eggsnew, gene_index::Int64, inputs::ExogenousInputs, t)

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)

    curr = get_lifestage(node, species, Egg)
    curr_μ, curr_q = temperature_effect(ctemp, curr)
    dens = compute_density(curr.density, E) #TODO: take curr

    prev = eggsnew

    dE[1,gene_index] = prev - E[1,gene_index]*(curr_μ * dens + curr_q * curr.n)
        for i in 2:size(dE)[1]
            dE[i,gene_index] = curr_q * curr.n * E[i-1,gene_index] - E[i,gene_index] * (curr_μ * dens + curr_q * curr.n)
        end

    return

end

"""
        function create_larva!(dL, L, E, node::Node, species::Type{<:Species}, gene_index::Int64, inputs, t)

    Returns larva state.
"""
function create_larva!(dL, L, E, node::Node, species::Type{<:Species}, gene_index::Int64, inputs::ExogenousInputs, t)

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)

    # TODO: UPDATE EACH TO "curr = stage_data"
    curr = get_lifestage(node, species, Larva)
    curr_μ, curr_q = temperature_effect(ctemp, curr)
    dens = compute_density(curr.density, L)

    prev = get_previous_lifestage(node, species, curr)
    prev_q = temperature_effect(ctemp, prev)[2]

    dL[1,gene_index] = prev_q * prev.n * E[end,gene_index] - L[1,gene_index] * (curr_μ * dens + curr_q * curr.n)

    for i in 2:size(dL)[1]
        dL[i,gene_index] = curr_q * curr.n * L[i-1,gene_index] - L[i,gene_index] * (curr_μ * dens + curr_q * curr.n)
    end

    # TODO: Change "dens" name to "compute_dens" or similar
    return

end

"""
        function create_pupa!(dP, P, L, node::Node, species::Type{<:Species}, gene_index::Int64, inputs, t)

    Returns pupa state.
"""
function create_pupa!(dP, P, L, node::Node, species::Type{<:Species}, gene_index::Int64, inputs::ExogenousInputs, t)

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)

    curr = get_lifestage(node, species, Pupa)
    curr_μ, curr_q = temperature_effect(ctemp, curr)

    dens = compute_density(curr.density, P)

    prev = get_previous_lifestage(node, species, curr)
    prev_q = temperature_effect(ctemp, prev)[2]

    dP[1,gene_index] = prev_q * prev.n * L[end, gene_index] - P[1,gene_index] * (curr_μ * dens + curr_q * curr.n)
    for i in 2:size(dP)[1]
        dP[i,gene_index] = curr_q * curr.n * P[i-1,gene_index] - P[i,gene_index] * (curr_μ * dens + curr_q * curr.n)
    end

    return

end

"""
        function create_male!(dM, M, P, Φ, Ξ_m, Ω, node::Node, species::Type{<:Species}, gene_index::Int64, inputs, t)

    Returns adult male state.
"""
function create_male!(dM, M, P, Φ, Ξ_m, Ω, node::Node, species::Type{<:Species}, gene_index::Int64, inputs::ExogenousInputs, t)

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)
    control = get_exogenous_intervention(inputs, node, species, Male, gene_index)

    curr = get_lifestage(node, species, Male)
    curr_μ, _ = temperature_effect(ctemp, curr)

    dens = compute_density(curr.density, M)

    prev = get_previous_lifestage(node, species, curr)
    prev_q = temperature_effect(ctemp, prev)[2]

    genetics = get_genetics(node, species)
    homozygous_modified = findfirst(isodd, genetics.all_modified)

    dM[gene_index] = ((1-Φ[gene_index]) * prev_q * prev.n * P[end,gene_index]) * Ξ_m[gene_index] - (1 + Ω[gene_index]) * curr_μ * M[gene_index] * dens + control

    return

end

"""
        function mate(P, M, Φ, Ξ_f, Η, node::Node, species::Type{<:Species}, gene_index::Int64, inputs, t)

    Mating.
"""
function mate(P, M, Φ, Ξ_f, Η, node::Node, species::Type{<:Species}, gene_index::Int64, inputs::ExogenousInputs, t)

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)

    prev = get_lifestage(node, species, Pupa)
    prev_q = temperature_effect(ctemp, prev)[2]

    nowmate = zeros(length(Η))
    for i in 1:length(Η)
        nowmate[i] = M[i]*Η[i]
    end
    nowmate = LinearAlgebra.normalize(nowmate, 1)

    mated = nowmate * (Φ[gene_index] * prev_q * prev.n * P[end,:]')

    return mated

end

"""
        function mate(P, M, Φ, Ξ_f, Η, node::Node, species::Type{<:Species}, gene_index::Int64, inputs, t)

    Returns adult female state.
"""
function create_female!(dF, F, Ω, Ξ_f, node::Node, species::Type{<:Species}, matematrix::Array, gene_index::Int64, inputs::ExogenousInputs, t)

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)

    curr = get_lifestage(node, species, Female)
    curr_μ, _ = temperature_effect(ctemp, curr)
    dens = compute_density(curr.density, F)

    for eachF in 1:size(dF)[1]
        control = get_exogenous_intervention(inputs, node, species, Female, gene_index, eachF)
        dF[eachF,gene_index] = matematrix[eachF, gene_index]*Ξ_f[eachF] - (1 + Ω[gene_index]) * curr_μ * F[eachF,gene_index] * dens + control
    end

    return

end