

function population_model_network(du, u, (network, inputs), t)

    ##################
    #   Parameters   #
    ##################

    n_nodes = count_nodes(network)

    for (index_node, node_pair) in enumerate(network.nodes)
        key_nodename, node = node_pair

        for (index_organism, key_species) in enumerate(keys(node.organisms))

            genetics = get_genetics(network, key_nodename, key_species)
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

            n_substages = count_substages(network, key_nodename, key_species)
            nE = n_substages[1]
            nL = n_substages[2]
            nP = n_substages[3]
            nJuv = nE+nL+nP
            nM = n_substages[4]
            nF = n_substages[5]*gN + nM + nJuv
    
            Q = get_migration(network, key_species)

            ##################
            #   State space  #
            ##################
            # TODO: Make indexing cleaner/more generalizable with: substages, state space (u, du), and Q.
            # Eg: Q could be a sparse matrix (all 0 in juv stages)

            E = u.x[index_node].x[index_organism][1:nE, :]
            dE = @view du.x[index_node].x[index_organism][1:nE, :]

            L = u.x[index_node].x[index_organism][1+nE : nE+nL, :]
            dL = @view du.x[index_node].x[index_organism][1+nE:nE+nL, :]

            P = u.x[index_node].x[index_organism][1+nE+nL : nE+nL+nP, :]
            dP = @view du.x[index_node].x[index_organism][1+nE+nL:nE+nL+nP, :]

            M = u.x[index_node].x[index_organism][1+nJuv, :]
            dM = @view du.x[index_node].x[index_organism][1+nJuv, :]

            F = u.x[index_node].x[index_organism][1+nM+nJuv : nF, :]
            dF = @view du.x[index_node].x[index_organism][1+nM+nJuv:nF, :]

            for gene_index in 1:gN

                eggsnew = oviposit(F, node, key_species, genetics, gene_index, inputs, t) 

                create_egg!(dE, E, node, key_species, eggsnew, gene_index, inputs, t)

                create_larva!(dL, L, E, node, key_species, gene_index, inputs, t)

                create_pupa!(dP, P, L, node, key_species, gene_index, inputs, t)

                create_male!(dM, M, P, Φ, Ξ_m, Ω, node, key_species, gene_index, inputs, t)

                matematrix = mate(P, M, Φ, Ξ_f, Η, node, key_species, gene_index, inputs, t)

                create_female!(dF, F, Ω, Ξ_f, node, key_species, matematrix, gene_index, inputs, t)

                ##################
                #    Migration   #
                ##################

                for n in 1:nE
                    dE[n, gene_index] = dE[n, gene_index] + Q[n, gene_index][index_node,:]'*[u.x[i].x[index_organism][n, gene_index] for i in 1:n_nodes]
                end

                for n in 1:nL
                    dL[n, gene_index] = dL[n, gene_index] + Q[n+nE, gene_index][index_node,:]'*[u.x[i].x[index_organism][n+nE, gene_index] for i in 1:n_nodes]
                end

                for n in 1:nP
                    dP[n, gene_index] = dP[n, gene_index] + Q[n+nE+nL, gene_index][index_node,:]'*[u.x[i].x[index_organism][n+nE+nL, gene_index] for i in 1:n_nodes]
                end

                dM[gene_index] = dM[gene_index] + Q[nM+nE+nL+nP, gene_index][index_node,:]'*[u.x[i].x[index_organism][nM+nE+nL+nP, gene_index] for i in 1:n_nodes]

                for n in 1:gN
                    dF[n, gene_index] = dF[n, gene_index] + Q[n+nE+nL+nP+nM, gene_index][index_node,:]'*[u.x[i].x[index_organism][n+nE+nL+nP+nM, gene_index] for i in 1:n_nodes]
                end

            end

        end

    end

end
