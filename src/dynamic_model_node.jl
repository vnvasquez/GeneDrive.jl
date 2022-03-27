
"""
    population_model_node(du, u, (network, inputs), t)

Solve node implementation of dynamic population model.
"""
function population_model_node(du, u, (node, inputs), t)

    ##################
    #   Parameters   #
    ##################

    for (index_organism, key_species) in enumerate(keys(node.organisms))
        genetics = get_genetics(node, key_species)
        gN = count_genotypes(genetics)
        cube = genetics.cube
        S = genetics.S
        Τ = genetics.Τ
        Φ = genetics.Φ
        Ξ_m = genetics.Ξ_m
        Ξ_f = genetics.Ξ_f
        Ω = genetics.Ω
        Β = genetics.Β
        Η = genetics.Η

        n = count_substages(node, key_species)
        nE = n[1]
        nL = n[2]
        nP = n[3]
        nJuv = nE + nL + nP
        nM = n[4]
        nF = n[5] * gN + nM + nJuv

        ##################
        #   State space  #
        # TODO: FLEXIBLY DEFINE STATE SPACE (AND SUBSTAGES) GIVEN ORGANISMS WITH VARIED BIOLOGY
        ##################

        E = u.x[index_organism][1:nE, :]
        dE = @view du.x[index_organism][1:nE, :]

        L = u.x[index_organism][(1 + nE):(nE + nL), :]
        dL = @view du.x[index_organism][(1 + nE):(nE + nL), :]

        P = u.x[index_organism][(1 + nE + nL):(nE + nL + nP), :]
        dP = @view du.x[index_organism][(1 + nE + nL):(nE + nL + nP), :]

        M = u.x[index_organism][1 + nJuv, :]
        dM = @view du.x[index_organism][1 + nJuv, :]

        F = u.x[index_organism][(1 + nM + nJuv):nF, :]
        dF = @view du.x[index_organism][(1 + nM + nJuv):nF, :]

        ##################
        #   Life stages  #
        ##################

        for gene_index in 1:gN
            eggsnew = oviposit(F, node, key_species, genetics, gene_index, inputs, t)

            create_egg!(dE, E, node, key_species, eggsnew, gene_index, inputs, t)

            create_larva!(dL, L, E, node, key_species, gene_index, inputs, t)

            create_pupa!(dP, P, L, node, key_species, gene_index, inputs, t)

            create_male!(dM, M, P, Φ, Ξ_m, Ω, node, key_species, gene_index, inputs, t)

            matematrix = mate(P, M, Φ, Ξ_f, Η, node, key_species, gene_index, inputs, t)

            create_female!(
                dF,
                F,
                Ω,
                Ξ_f,
                node,
                key_species,
                matematrix,
                gene_index,
                inputs,
                t,
            )
        end
    end
end
