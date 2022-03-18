
# include("paper_enviro_data.jl")
using JuMP, Ipopt
using RecursiveArrayTools
using NLsolve
using DataStructures
#using KNITRO
include("src/GeneDrive.jl")
# TODO: to facilitate laptop work with Pardiso, see paper_optimization_initialexploration.jl in /TestDrive

##################
#      Node      #
##################

organisms = OrderedDict(AedesAegypti => Organism{AedesAegypti}(
    #genetics_wolbachia(),
    genetics_ridl(),
    #stages_noresponse_500())); TODO: check this
    stages_rossi_500()));

testlab = fill(27.0, 365);
node = Node(:Cuenca, organisms, TimeSeriesTemperature(testlab), (2.9001, 79.0059));

initial_condition, density_parameter = init_node!(node);
initial_condition.x[1]
density_parameter[1][1] # for logistic use Kappa

##################
#      Data      #
##################

organism_count = count_organisms(node)
gene_count = count_genotypes(node, AedesAegypti)
stage_count = sum(count_substages(node, AedesAegypti)[1:4]) + count_substages(node, AedesAegypti)[5]*gene_count

wildtype = get_wildtype(node, AedesAegypti)
homozygous_modified = get_homozygous_modified(node, AedesAegypti)

ctemp = get_initial_temperature(node)
species = AedesAegypti

# Genetics
cube = node.organisms[species].gene_data.cube
Β = node.organisms[species].gene_data.Β
Η = node.organisms[species].gene_data.Η
Τ = node.organisms[species].gene_data.Τ
Φ = node.organisms[species].gene_data.Φ
S = node.organisms[species].gene_data.S
Ω = node.organisms[species].gene_data.Ω
Ξ_m = node.organisms[species].gene_data.Ξ_m
Ξ_f = node.organisms[species].gene_data.Ξ_f

# Lifestages w/ Temperature
egg = get_lifestage(node, species, Egg)
μE, qE = temperature_effect(ctemp, egg)
nE = egg.n
densE = egg.density

larva = get_lifestage(node, species, Larva)
μL, qL = temperature_effect(ctemp, larva)
nL = larva.n
densL = larva.density

pupa = get_lifestage(node, species, Pupa)
μP, qP = temperature_effect(ctemp, pupa)
nP = pupa.n
densP = pupa.density

male = get_lifestage(node, species, Male)
μM, _ = temperature_effect(ctemp, male)
nM = male.n
densM = male.density

female = get_lifestage(node, species, Female)
μF, _ = temperature_effect(ctemp, female)
nF = female.n*gene_count
densF = female.density

tspan = (1.0, 365.0)
tempseries = TemperatureSeriesData(node, collect(tspan[1]:tspan[2]), testlab);

##################
#      Sets      #
##################

# Time
T = 1:2

# Genetics
G = 1:gene_count

# Stages/substages
SE = 1:nE
SL = 1:nL
SP = 1:nP
SM = 1:nM
SF = 1:nF

# Initial conditions mapped to stages TODO: check this
#=
SE_map = SE
SL_map = nE+1 : nE+nL
SP_map = nE+nL+1 : nE+nL+nP
SM_map = nE+nL+nP+1
SF_map = nE+nL+nP+nM+1 : nE+nL+nP+nM+nF
=#

# Organisms
O = 1:organism_count

##################
#    Solver      # Linear_solver:  HSL/ma86 = iMac || Pardiso = laptop
##################

i = optimizer_with_attributes(Ipopt.Optimizer, "linear_solver" => "Pardiso") # "linear_solver" => "Pardiso" #"linear_solver" => "ma86"
#i = optimizer_with_attributes(KNITRO.Optimizer, "par_lsnumthreads" => 8, "honorbnds" => 1, "linsolver"=> 8, "algorithm" => 1) TODO: check this

model = Model(i);

###########################################
# DECLARE VARIABLES_1: Life Stages
###########################################

@variable(model, E[O, SE, G, T] >= 0)
@variable(model, L[O, SL, G, T] >= 0)
@variable(model, P[O, SP, G, T] >= 0)
@variable(model, M[O, SM, G, T] >= 0)
@variable(model, F[O, SF, G, T] >= 0)

###########################################
# WARMSTART VARIABLES_1
###########################################

# Access initial conditions
E0 = initial_condition.x[1][SE,:]
L0 = initial_condition.x[1][nE+1:nL+nE,:]
P0 = initial_condition.x[1][nE+nL+1:nE+nL+nP,:]
M0 = initial_condition.x[1][nE+nL+nP+1,:]' # NB transpose
F0 = initial_condition.x[1][nE+nL+nP+2:end,:]

# Warmstart (option 2): Furnish start values for *EVERY* timestep to reduce search time
for t in T
    set_start_value.(E[1, :,:,t].data, E0)
    set_start_value.(L[1, :,:,t].data, L0)
    set_start_value.(P[1, :,:,t].data, P0)
    set_start_value.(M[1, :,:,t].data, M0)
    set_start_value.(F[1, :,:,t].data, F0)
end

##########################################
# DECLARE VARIABLES_2: Restrict and otherwise schedule controls
# also see CONSTRAINTS_B (OPERATIONAL CONSTRAINTS: Sets overall limit on releases)
###########################################

# Bound the control variable (size of release permitted per time step) TODO: currently applicable only to adult stages; update for PGSIT (eggs)
@variable(model, 0.0 <= control_M[O, SM, G, T]) # <= 5000.0) # REMOVE "upper" part of the bound when using a schedule/envelope
@variable(model, 0.0 <= control_F[O, SF, G, T]) # <= 5000.0)

# Clear: Restrict releases such that none allowed for any organisms/stages/genotypes/timesteps
fix.(control_M, 0.0; force = true)
fix.(control_F, 0.0; force = true) # activate when checking dynamics
fix.(control_F[:, wildtype, :, :], 0.0; force = true) # use always bc F = matrix
fix.(control_F[:, :, homozygous_modified, :], 0.0; force = true) # use always bc F = matrix

### This sets up releases to check dynamics
# Important: Specify in which nodes/orgs/genes/timesteps the releases are permitted
# FEMALES
#=
for t in T[25:10:65]
        @show t
        fix.(control_F[1, SF[homozygous_modified], G[wildtype], t], 50.0; force = true)
end
=#
# MALES
#=
for t in T[25:10:65]
    @show t
    fix.(control_M[1, :, G[homozygous_modified], t], 50.0; force = true)
end
=#


###########################################
# CONSTRAINTS_A: BIOLOGICAL CONSTRAINTS
###########################################

#### EGGS
# First state space, first time step, all genes
@constraint(model, E_con_A0[o in O, s = SE[1], g = G, t = T[1]],

        # Initial condition (time = 0 so not indexed here)
        E[o, s, g, t] == E0[s, g] +

        # Look back to oviposition (aka last substage of previous stage, current timestep)
        sum((cube[:,:,g].*Τ[:,:,g].*S[g]*Β[g].*F[o,:,:,t].data)[:,:]) -

        # Change: die or move to next substage
        E[o, s, g, t] * (μE*compute_density(densE, sum(E[o,:,:,t])) + qE*nE))

# First state space, second through last time step, all genes
@constraint(model, E_con_A1[o in O, s = SE[1], g = G, t = T[2:end]],

        # First state space, look back to previous time step (time = 1)
        E[o, s, g, t] ==  E[o, s, g, t-1] +

        # Look back to oviposition (aka last substage of previous stage, current timestep)
        sum((cube[:,:,g].*Τ[:,:,g].*S[g]*Β[g].*F[o,:,:,t].data)[:,:]) -

        # Change: die or move to next substage
        E[o, s, g, t] * (μE*compute_density(densE, sum(E[o,:,:,t])) + qE*nE))

# Second through last state space, first time step, all genes
@constraint(model, E_con_B0[o in O, s = SE[2:end], g = G, t = T[1]],

        # Initial condition (time = 0 so not indexed here) TODO: CLARIFY THIS
        E[o, s, g, t] == E0[ s, g] +

        # Look back to previous state space, current timestep
        qE*nE*E[o, s-1, g, t] -

        # Change: die or move to next substage
        E[o, s, g, t] * (μE*compute_density(densE, sum(E[o,:,:,t])) + qE*nE))

# Second through last state space, second through last time step, all genes
@constraint(model, E_con_B1[o in O, s = SE[2:end], g = G, t = T[2:end]] ,

        # Start at second state space, look back to previous time step
        E[o, s, g, t] ==   E[o, s, g, t-1] +

        # Look back to previous state space, current timestep
        qE*nE*E[o, s-1, g, t] -

        # Change: die or move to next STAGE
        E[o, s, g, t] * (μE*compute_density(densE, sum(E[o,:,:,t])) + qE*nE))


#### LARVAE
# First state space, first time step, all genes
@constraint(model, L_con_A0[o in O, s = SL[1], g = G, t = T[1]],

        # Initial condition (time = 0 so not indexed here)
        L[o, s, g, t] ==  L0[ s, g] +

        # Look back to last substage of previous stage, current timestep
        qE * nE * E[o,end,g,t] -

        # Change: die or move to next substage
        L[o, s, g, t] * (μL*compute_density(densL, sum(L[o,:,:,t])) + qL*nL))

# First state space, second through last timestep, all genes
@constraint(model, L_con_A1[o in O, s = SL[1], g = G, t = T[2:end]],

        # Start at first substage, look back to previous timestep
         L[o, s, g, t] ==  L[o, s, g, t-1] +

        # Look back to the last substage of previous stage, current timestep
        qE * nE * E[o, end,g,t] -

        # Change: die or move to the next substage
        L[o, s, g, t] * (μL*compute_density(densL, sum(L[o,:,:,t])) + qL*nL))

# Second through last state space, first time step, all genes
@constraint(model, L_con_B0[o in O, s = SL[2:end], g = G, t = T[1]],

        # Initial condition (time = 0 so not indexed here)
        L[o, s, g, t] == L0[ s, g] +

        # Look back to the previous substage, current timestep
        qL*nL*L[o, s-1, g, t] -

        # Change: die or move to the next substage
        L[o, s, g, t] * (μL*compute_density(densL, sum(L[o,:,:,t])) + qL*nL))

# Second through last state space, second through last time step, all genes
@constraint(model, L_con_B1[o in O, s = SL[2:end], g = G, t = T[2:end]],

        # Start at second substage, look back to previous timestep
        L[o, s, g, t] == L[o, s, g, t-1] +

        # Look back to the previous substage, current timestep
         qL*nL*L[o, s-1, g, t] -

        # Change: die or move to the next STAGE
         L[o, s, g, t] * (μL*compute_density(densL, sum(L[o,:,:,t])) + qL*nL))


#### PUPAE
# First state space, first time step, all genes
@constraint(model, P_con_A0[o in O, s = SP[1], g = G, t = T[1]],

        # Initial condition (time = 0 so not indexed here)
        P[o, s, g, t] ==  P0[ s, g] +

        # Look back to the last substage of the previous stage
        qL * nL * L[o, end , g ,t] -

        # Change: die or move to the next substage
        P[o, s, g, t] * (μP*compute_density(densP, sum(P[o, :, :, t])) + qP*nP))

# First state space, second through last time step, all genes
@constraint(model, P_con_A1[o in O, s = SP[1], g = G, t = T[2:end]],

        # Start at first substage, look back to previous timestep
        P[o, s, g, t] ==  P[o, s, g, t-1] +

        # Look back to the last substage of the previous stage
        qL * nL * L[o, end , g ,t] -

        # Change: die or move to the next substage
        P[o, s, g, t] * (μP*compute_density(densP, sum(P[o, :, :, t])) + qP*nP))

# Second through last state space, first time step, all genes
@constraint(model, P_con_B0[o in O, s = SP[2:end], g = G, t = T[1]] ,

        # Initial condition (time = 0 so not indexed here)
        P[o, s, g, t] ==  P0[ s, g] +

        # Look back to previous substage, current timestep
        qP*nP*P[o, s-1, g, t] -

        # Change: die or move to the next substage
        P[o, s, g, t] * (μP*compute_density(densP, sum(P[o, :, :, t])) + qP*nP))

# Second through last state space, second through last time step, all genes
@constraint(model, P_con_B1[o in O, s = SP[2:end], g = G, t = T[2:end]],

        # Start at second substage, look back to the previous timestep
        P[o, s, g, t] ==  P[o, s, g, t-1] +

        # Look back to the previous substage, current time step
        qP*nP*P[o, s-1, g, t] -

        # Change: die or move to the next STAGE
        P[o, s, g, t] * (μP*compute_density(densP, sum(P[o, :, :, t])) + qP*nP))

#### MALES
# First state space, first timestep
@constraint(model, M_con_0[o in O, s = SM, g = G, t = T[1]],

        # Initial condition (time = 0 so not indexed here)
        M[o, s, g, t] ==  M0[ s, g] +

        # Look back to last substage pf previous STAGE and take sex-appropriate (phi) portion
        (1-Φ[g]) * qP * nP * P[o, end, g ,t] * Ξ_m[g] -

        # Change: die or move to next STAGE
        (1 + Ω[g]) * μM * M[o, s, g, t] * compute_density(densM, sum(M[o, :, :, t])))

# First state space, second through last timestep
@constraint(model, M_con_1[o in O, s = SM, g = G, t = T[2:end]],

        # First substage, look back to previous timstep
        M[o, s, g, t] ==  M[o, s, g, t-1] +

        # Look back to last substage pf previous STAGE and take sex-appropriate (phi) portion
        (1-Φ[g]) * qP * nP * P[o, end, g ,t] * Ξ_m[g] -

        # Change: die or move to next STAGE + CONTROL when applicable
        (1 + Ω[g]) * μM * M[o, s, g, t] * compute_density(densM, sum(M[o, :, :, t])) + control_M[o, s, g, t])


# Mating constraint to save computation time
@constraint(model, mate_bound[o in O, g = G, t = T], M[o, 1, g, t]*Η[g] <= (sum(M[o, 1, i, t]*Η[i] for i in G)))


#### FEMALES
# First state space, first timestep
@NLconstraint(model, F_con_0[o in O, s in SF, g in G, t = T[1]],

        # Initial condition (time = 0 so not indexed here)
        F[o, s, g, t] == F0[ s, g] +

        # Mating term (the only nonlinearity) multipled by sex-appropriate (phi) portion of previous STAGE
        (M[o, 1, g, t]*Η[g]/(1e-6 + sum(M[o, 1, i, t]*Η[i] for i in G)))*(Φ[s] * qP * nP * P[o, end, s, t] * Ξ_f[s]) -

        # Change: die or move to next STAGE (oviposition/egg production)
        (1 + Ω[g]) * μF * F[o, s, g, t] )

# First state space, second through last timesteps
@NLconstraint(model, F_con_1[o in O, s in SF, g in G, t in T[2:end]],

        # First stage, look back to previous timestep
        F[o, s, g, t] == F[o, s, g, t-1] +

        # Mating term (the only nonlinearity) multipled by sex-appropriate (phi) portion of previous STAGE
        (M[o, 1, g, t]*Η[g]/(1e-6+sum(M[o, 1, i, t]*Η[i] for i in G)))*(Φ[s] * qP * nP * P[o, end, s, t] * Ξ_f[s]) -

        # Change: die or move to next STAGE (oviposition/egg production) + CONTROL when applicable
        (1 + Ω[g]) * μF * F[o, s, g, t] + control_F[o, s, g, t])


###########################################
# CONSTRAINTS_B: OPERATIONAL CONSTRAINTS
###########################################


#TODO: COME BACK TO THESE


###########################################
# OBJECTIVE FUNCTIONS: Use one at a time
###########################################

#### OPTION 1 [DEBUGGING]
@objective(model, Min, 0)

#TODO: COME BACK TO THESE


###########################################
# OPTIMIZE!
###########################################

# Optimize using given objective function (can also use without ANY specific objective function)
optimize!(model)

# Nice way to monitor what's going on
println("Objective value, single node multi-organism: ", objective_value(model))
