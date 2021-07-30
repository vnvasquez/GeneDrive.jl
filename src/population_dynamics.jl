################################################################################
#                         Population: ODE Formulation                          #
################################################################################

"""
        function oviposit(F, cube, Τ, S, Β, gene_index::Int64)

    Returns oviposited eggs (count).
"""
function oviposit(F, node::Node, key_species, genetics::Genetics{C}, gene_index::Int64, inputs::ExogenousInputs, t) where {C <: Construct}
#function oviposit(F, node::Node, genetics, gene_index::Int64, t)

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)

    cube = genetics.cube
    Τ = genetics.Τ
    S = genetics.S
    Β = genetics.Β

    #TODO: This needs updating since it allocates a matrix per iteration
    ΒS = Matrix{Float64}(undef, length(Β), length(S))
    for i in 1:length(Β)
        ΒS[i,:] = Β[i].*S'
    end
    # TODO: add genetic parameters necessary for additional constructs
    O = cube[:,:,gene_index].*Τ[:,:,gene_index].*ΒS.*F

    return sum(O)
end

"""
        function get_response_type(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})

    Returns temperature response type for appropriate dispatch of temperature-sensitive functions.
"""
function get_response_type(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})
    return typeof(node.organisms[species].all_stages[life_stage].μ_temperature_response)
end

"""
        function temperature_responsive_inheritence(response_type::Type{NoResponse}, ctemp::Float64, genetics::Genetics, gene_index::Int64)

    Dispatches on `NoResponse` type to return offspring likelihoods for each genotype that are not responsive to temperature.
"""
function temperature_responsive_inheritence(response_type::Type{NoResponse}, ctemp::Float64, genetics::Genetics, gene_index::Int64)
    return genetics.cube[:,:,gene_index]
end

"""
        function temperature_responsive_inheritence(response_type::Type{T}, ctemp::Float64, genetics::Genetics, gene_index::Int64) where {T <:TemperatureResponse}

    Dispatches on `TemperatureResponse` type(s) to return offspring likelihoods for each genotype that are responsive to temperature.
"""
function temperature_responsive_inheritence(response_type::Type{T}, ctemp::Float64, genetics::Genetics, gene_index::Int64) where {T <:TemperatureResponse}
    # TODO: citation: Fig 7B + associated data for Ross et al (2019)
    if ctemp >= 35.0 && gene_index == 1
        return genetics.cube[:,:,gene_index] = [0.0 0.0; 0.0 0.0]
    elseif ctemp >= 35.0 && gene_index == 2
        return genetics.cube[:,:,gene_index] = [1.0 1.0; 1.0 1.0]
    else
        return genetics.cube[:,:,gene_index]
    end
end

"""
        function temperature_responsive_hatching(response_type::Type{NoResponse}, ctemp::Float64, genetics::Genetics, gene_index::Int64)

    Dispatches on `NoResponse` type to return hatch proportions for each genotype that are not responsive to temperature.
"""
function temperature_responsive_hatching(response_type::Type{NoResponse}, ctemp::Float64, genetics::Genetics, gene_index::Int64)
    return genetics.Τ[:,:,gene_index]
end

"""
        function temperature_responsive_hatching(response_type::Type{NoResponse}, ctemp::Float64, genetics::Genetics, gene_index::Int64)

    Dispatches on `TemperatureResponse` type(s) to return hatch proportions for each genotype that are responsive to temperature.
"""
function temperature_responsive_hatching(response_type::Type{T}, ctemp::Float64, genetics::Genetics, gene_index::Int64) where {T <:TemperatureResponse}
    # TODO: citation: Fig 7A + associated data for Ross et al (2019)
    if gene_index == 2
        return genetics.Τ[:,:,gene_index]
    else
        hatch = 1.0
        if ctemp <= 26.0
            hatch = 0.91833
        elseif ctemp > 36.0
            hatch = 0.0
        else
            hatch = 0.91832971025 - 0.021*(ctemp - 26) - 0.0055*exp(0.485*(ctemp - 26))
        end
        updated_tau = genetics.Τ[:,:,gene_index]
        updated_tau[1,:] .= hatch
    end

    # TODO: verify whether CI should be occuring in both gene_index 1 & 2 ->
    # seems only applicable to gene_1 (but are C & G entries - and therefore K & J - even correct?
    if ctemp >= 35.0 && gene_index == 1
        return updated_tau[2,1] = 1.0
    elseif ctemp >= 35.0 && gene_index == 2
        return updated_tau[2,1] = 1.0 #TODO: check to see if it is necessary to differentiate tau_1 vs tau_2; if not collabse into single if/else without "&& gene_index" piece
    else
        return updated_tau
    end
end

"""
        function oviposit(F, node::Node, genetics::Genetics{Wolbachia}, gene_index::Int64, inputs::ExogenousInputs, t)

    Returns oviposited eggs (count). Includes capacity for temperature-sensitive hatching and offspring likelihoods.
"""
function oviposit(F, node::Node, key_species, genetics::Genetics{Wolbachia}, gene_index::Int64, inputs::ExogenousInputs, t) 

    node_name = get_name(node)
    ctemp = get_temperature_value(node.temperature, inputs.temperature[node_name], t)
    response_type = get_response_type(node, key_species, Egg)

    cube = temperature_responsive_inheritence(response_type, ctemp, genetics, gene_index)
    Τ = temperature_responsive_hatching(response_type, ctemp, genetics, gene_index)
    S = genetics.S
    Β = genetics.Β

    ΒS = Matrix{Float64}(undef, length(Β), length(S))
    for i in 1:length(Β)
        ΒS[i,:] = Β[i].*S'
    end

    O = cube.*Τ.*ΒS.*F
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

    @assert !all(isnan.(nowmate)) Η, M

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
        dF[eachF, gene_index] = matematrix[gene_index, eachF]*Ξ_f[eachF] - (1 + Ω[eachF]) * curr_μ * F[eachF,gene_index] * dens + control # two indexing corrections
    end

    return

end
