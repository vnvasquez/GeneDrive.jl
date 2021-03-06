################################################################################
#             Species and Stage-Specific Temperature Responses                 #
################################################################################

abstract type TemperatureResponse end

"""
    get_temperature_response(ctemp::Float64, response::TemperatureResponse)

Return `TemperatureResponse` data (model with temperature response).
"""
function get_temperature_response(ctemp::Float64, response::TemperatureResponse)
    return get_temperature_response(ctemp::Float64, response::TemperatureResponse, 0.0)
end

"""
    mutable struct NoResponse <: TemperatureResponse
        baseline_value::Float64
    end

Data for model without temperature response.

# Arguments:

  - `baseline_value::Float64`: Literature-sourced (rather than dynamically calculated) vital rate parameter.
"""
mutable struct NoResponse <: TemperatureResponse
    baseline_value::Float64
end

"""
    get_temperature_response(::Float64, response::NoResponse, ::Float64) = response.baseline_value

Return `NoResponse` data (model without temperature response).
"""
get_temperature_response(::Float64, response::NoResponse, ::Float64) =
    response.baseline_value

"""
    temperature_effect(ctemp::Float64, stage)

Return effect of temperature on organism vital rates (mortality, duration = `μ_temperature_response`, `q_temperature_response`).
"""
function temperature_effect(ctemp::Float64, stage)
    μ_temperature_response = stage.μ_temperature_response
    q_temperature_response = stage.q_temperature_response
    duration = get_temperature_response(ctemp, q_temperature_response)
    mortality = get_temperature_response(ctemp, μ_temperature_response, duration)
    return mortality, duration
end

#TODO: include validation range for each lit-based parameterization (validation range: `(0, nothing)`).
# TODO: Complete docstrings with inclusion of (1) baseline struct and (2) argument description

########################################
#                Rossi                 #
########################################

"""
    mutable struct EggMortalityRossi <: TemperatureResponse
        a::Float64
        b::Float64
    end

Data for `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Egg. Source: Rossi et al (2014).

# Arguments

  - `a::Float64`: Death rate at low temperatures, validation range: `(0, nothing)`
  - `b::Float64`: Influence factor of temperature, validation range: `(0, nothing)`
"""
mutable struct EggMortalityRossi <: TemperatureResponse
    a::Float64
    b::Float64
end

"""
    get_temperature_response(ctemp::Float64, response::EggMortalityRossi, ::Float64)

Return `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Egg. Source: Rossi et al (2014).
"""
function get_temperature_response(ctemp::Float64, response::EggMortalityRossi, ::Float64)
    return response.a * exp(response.b * ctemp)
end

"""
    mutable struct EggDurationRossi <: TemperatureResponse
        a::Float64
        b::Float64
        c::Float64
        d::Float64
        e::Float64
        f::Float64
    end

Data for `q_temperature_response`. `Species`: AedesAegypti, `LifeStage`:, egg stage. Source: Rossi et al (2014).
"""
mutable struct EggDurationRossi <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
end

"""
    get_temperature_response(ctemp::Float64, response::EggDurationRossi, ::Float64)

Return `q_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Egg. Source: Rossi et al (2014).
"""
function get_temperature_response(ctemp::Float64, response::EggDurationRossi, ::Float64)
    return response.a *
           (ctemp + response.b) *
           (
               exp(response.c - (response.d / (ctemp + response.b))) / 1 +
               exp(response.e - (response.f / (ctemp + response.b)))
           )
end

"""
Data for `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Larva. Source: Rossi et al (2014).
"""
mutable struct LarvaMortalityRossi <: TemperatureResponse
    a::Float64
    b::Float64
    # TODO: why temp dependent deathrate unused in source?
    # c:: Float64 # π0 = 4.03*10^(-6)
end

"""
    get_temperature_response(ctemp::Float64, response::LarvaMortalityRossi, ::Float64)

Return `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`:, larval stage. Source: Rossi et al (2014).
"""
function get_temperature_response(ctemp::Float64, response::LarvaMortalityRossi, ::Float64)
    return response.a * exp(response.b * ctemp^2)
end

"""
Data for `q_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Larva. Source: Rossi et al (2014).
"""
mutable struct LarvaDurationRossi <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
end

"""
    get_temperature_response(ctemp::Float64, response::LarvaDurationRossi, ::Float64)

Return `q_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Larva. Source: Rossi et al (2014).
"""
function get_temperature_response(ctemp::Float64, response::LarvaDurationRossi, ::Float64)
    return response.a *
           (ctemp + response.b) *
           exp(response.c - (response.d / (ctemp + response.b)))
end

"""
Data for `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Pupa. Source: Rossi et al (2014).
"""
mutable struct PupaMortalityRossi <: TemperatureResponse
    a::Float64
    b::Float64
    # TODO: why temp dependent deathrate unused in source?
    # c:: Float64 # π0 = 4.03*10^(-6)
end

"""
    get_temperature_response(ctemp::Float64, response::PupaMortalityRossi, ::Float64)

Return `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Pupa. Source: Rossi et al (2014).
"""
function get_temperature_response(ctemp::Float64, response::PupaMortalityRossi, ::Float64)
    return response.a * exp(response.b * ctemp^2)
end

"""
Data for `q_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Pupa. Source: Rossi et al (2014) and Poletti et al (2011) Table 1.
"""
mutable struct PupaDurationRossi <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
end

"""
    get_temperature_response(ctemp::Float64, response::PupaDurationRossi, ::Float64)

Return `q_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Pupa. Source: Rossi et al (2014) and Poletti et al (2011) Table 1.
"""
function get_temperature_response(ctemp::Float64, response::PupaDurationRossi, ::Float64)
    return response.a * ctemp^2 + response.b * ctemp + response.c
end

"""
Data for `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Male, Female. Source: Rossi et al (2014).
"""
mutable struct AdultMortalityRossi <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
end

"""
    get_temperature_response(ctemp::Float64, response::AdultMortalityRossi, ::Float64)

Return `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Male, Female. Source: Rossi et al (2014).
"""
function get_temperature_response(ctemp::Float64, response::AdultMortalityRossi, ::Float64)
    return response.a * exp(-response.b * (ctemp - response.c)^2) +
           response.d * (ctemp - response.c)^2
end

########################################
#             El Moustaid              #
########################################
"""
Data for `q_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Egg. Source: El Moustaid et al (2019).
"""
mutable struct EggDurationMoustaid <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
end

"""
    get_temperature_response(ctemp::Float64, response::EggDurationMoustaid, ::Float64)

Return `q_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Egg. Source: El Moustaid et al (2019).
"""
function get_temperature_response(ctemp::Float64, response::EggDurationMoustaid, ::Float64)
    # TODO: why q = 1/duration in original implementation (applies to all juv)?
    return 1 / (response.a - (response.b * ctemp) + response.c * ctemp^2)
end

"""
Data for `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Egg. Source: El Moustaid et al (2019).
"""
mutable struct EggMortalityMoustaid <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
end

"""
    get_temperature_response(ctemp::Float64, response::EggMortalityMoustaid, duration::Float64)

Return `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Egg. Source: El Moustaid et al (2019).
"""
function get_temperature_response(
    ctemp::Float64,
    response::EggMortalityMoustaid,
    duration::Float64,
)
    survival = 0.0
    if ctemp > response.a || ctemp < response.b
        survival = 0.0
    elseif ctemp > response.a && ctemp <= response.a
        survival = -ctemp * response.c + response.d
    else
        survival = response.e * ctemp * (ctemp - response.b) * ((response.a - ctemp)^0.5)
    end
    return (1 - survival) / (1 / duration)
end

"""
Data for `q_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Larva. Source: El Moustaid et al (2019).
"""
mutable struct LarvaDurationMoustaid <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
end

"""
    get_temperature_response(ctemp::Float64, response::LarvaDurationMoustaid, ::Float64)

Return `q_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Larva. Source: El Moustaid et al (2019).
"""
function get_temperature_response(
    ctemp::Float64,
    response::LarvaDurationMoustaid,
    ::Float64,
)
    return 1 / (response.a - (response.b * ctemp) + response.c * ctemp^2)
end

"""
Data for `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Larva. Source: El Moustaid et al (2019).
"""
mutable struct LarvaMortalityMoustaid <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
end

"""
    get_temperature_response(ctemp::Float64, response::LarvaMortalityMoustaid, duration::Float64)

Return `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Larva. Source: El Moustaid et al (2019).
"""
function get_temperature_response(
    ctemp::Float64,
    response::LarvaMortalityMoustaid,
    duration::Float64,
)
    survival = 0.0
    if ctemp > response.a || ctemp < response.b
        survival = 0.0
    elseif ctemp > response.a && ctemp <= response.a
        survival = -ctemp * response.c + response.d
    else
        survival = response.e * ctemp * (ctemp - response.b) * ((response.a - ctemp)^0.5)
    end
    return (1 - survival) / (1 / duration)
end

"""
Data for `q_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Pupa. Source: El Moustaid et al (2019).
"""
mutable struct PupaDurationMoustaid <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
end

"""
    get_temperature_response(ctemp::Float64, response::PupaDurationMoustaid, ::Float64)

Return `q_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Pupa. Source: El Moustaid et al (2019).
"""
function get_temperature_response(ctemp::Float64, response::PupaDurationMoustaid, ::Float64)
    return 1 / (response.a - (response.b * ctemp) + response.c * ctemp^2)
end

"""
Data for `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Pupa. Source: El Moustaid et al (2019).
"""
mutable struct PupaMortalityMoustaid <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64 # El Moustaid personal communication: lit value incorrect
    d::Float64
    e::Float64
end

"""
    get_temperature_response(ctemp::Float64, response::PupaMortalityMoustaid, duration::Float64)

Return `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Pupa. Source: El Moustaid et al (2019).
"""
function get_temperature_response(
    ctemp::Float64,
    response::PupaMortalityMoustaid,
    duration::Float64,
)
    survival = 0.0
    if ctemp > response.a || ctemp < response.b
        survival = 0.0
    elseif ctemp > response.a && ctemp <= response.a
        survival = -ctemp * response.c + response.d
    else
        survival =
            min(1, response.e * ctemp * (ctemp - response.b) * ((response.a - ctemp)^0.5))
    end
    return (1 - survival) / (1 / duration)
end

"""
Data for `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Male, Female. Source: El Moustaid et al (2019).
"""
mutable struct AdultMortalityMoustaid <: TemperatureResponse
    a::Float64 # El Moustaid personal comm: lit value incorrect
    b::Float64 # El Moustaid personal comm: lit value incorrect
    c::Float64 # El Moustaid personal comm: use "lf" param from Mordecai 2017, Table B
end

"""
    get_temperature_response(ctemp::Float64, response::AdultMortalityMoustaid, ::Float64)

Return `μ_temperature_response`. `Species`: AedesAegypti, `LifeStage`: Male, Female. Source: El Moustaid et al (2019).
"""
function get_temperature_response(
    ctemp::Float64,
    response::AdultMortalityMoustaid,
    ::Float64,
)
    survival = 0.0
    if ctemp > response.a || ctemp < response.b
        survival = 0.0
    else
        survival = -response.c * (ctemp - response.b) * (response.a - ctemp)
    end

    μ = 0.0
    if ctemp > response.a || ctemp < response.b
        μ = 1.0
    else
        μ = 1 / survival
    end

    return μ
end

########################################
#               Abiodun                #
########################################
"""
Data for `q_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Egg. Source: Abiodun et al (2016).
"""
mutable struct EggDurationAbiodun <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
end

"""
Data for `μ_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Egg. Source: Abiodun et al (2016).
"""
mutable struct EggMortalityAbiodun <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
end

# Mortality and duration use the same coefficients
function _get_coefficients(
    ctemp::Float64,
    response::Union{EggDurationAbiodun, EggMortalityAbiodun},
)
    ctemp += response.a
    q =
        response.b +
        response.c * (response.d + (ctemp / response.e)^response.f)^(-response.g)
    return q
end

"""
    get_temperature_response(ctemp::Float64, response::EggDurationAbiodun, ::Float64)

Return `q_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Egg. Source: Abiodun et al (2016).
"""
function get_temperature_response(ctemp::Float64, response::EggDurationAbiodun, ::Float64)
    return _get_coefficients(ctemp, response)
end

"""
    get_temperature_response(ctemp::Float64, response::EggMortalityAbiodun, ::Float64)

Return `q_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Egg. Source: Abiodun et al (2016).
"""
function get_temperature_response(ctemp::Float64, response::EggMortalityAbiodun, ::Float64)
    q = _get_coefficients(ctemp, response)
    return exp(-1 / q)
end

"""
Data for `q_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Larva. Source: Abiodun et al (2016).
"""
mutable struct LarvaDurationAbiodun <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
    h::Float64
    i::Float64
end

"""
Data for `μ_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Larva. Source: Abiodun et al (2016).
"""
mutable struct LarvaMortalityAbiodun <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
    h::Float64
    i::Float64
end

# Mortality and duration use the same coefficients
function _get_coefficients(
    ctemp::Float64,
    response::Union{LarvaDurationAbiodun, LarvaMortalityAbiodun},
)
    ctemp += response.a
    qE =
        response.h +
        response.i * (response.d + (ctemp / response.e)^response.f)^(-response.g)
    qL =
        response.b +
        response.c * (response.d + (ctemp / response.e)^response.f)^(-response.g)
    return qE, qL
end

"""
    get_temperature_response(ctemp::Float64, response::LarvaDurationAbiodun, ::Float64)

Return `q_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Larva. Source: Abiodun et al (2016).
"""
function get_temperature_response(ctemp::Float64, response::LarvaDurationAbiodun, ::Float64)
    qE, qL = _get_coefficients(ctemp, response)
    return qL - qE
end

"""
    get_temperature_response(ctemp::Float64, response::LarvaMortalityAbiodun, ::Float64)

Return `q_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Larva. Source: Abiodun et al (2016).
"""
function get_temperature_response(
    ctemp::Float64,
    response::LarvaMortalityAbiodun,
    ::Float64,
)
    _, qL = _get_coefficients(ctemp, response)
    return exp(-1 / qL)
end

"""
Data for `q_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Pupa. Source: Abiodun et al (2016).
"""
mutable struct PupaDurationAbiodun <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
    h::Float64
    i::Float64
    j::Float64
    k::Float64
end

"""
Data for `μ_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Pupa. Source: Abiodun et al (2016).
"""
mutable struct PupaMortalityAbiodun <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
    h::Float64
    i::Float64
    j::Float64
    k::Float64
end

# Mortality and duration use the same coefficients
function _get_coefficients(
    ctemp::Float64,
    response::Union{PupaDurationAbiodun, PupaMortalityAbiodun},
)
    ctemp += response.a
    qL =
        response.h +
        response.i * (response.d + (ctemp / response.j)^response.k)^(-response.g)
    qP =
        response.b +
        response.c * (response.d + (ctemp / response.e)^response.f)^(-response.g)
    return qL, qP
end

"""
    get_temperature_response(ctemp::Float64, response::PupaDurationAbiodun, ::Float64)

Return `q_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Pupa. Source: Abiodun et al (2016).
"""
function get_temperature_response(ctemp::Float64, response::PupaDurationAbiodun, ::Float64)
    qL, qP = _get_coefficients(ctemp, response)
    return qP - qL
end

"""
    get_temperature_response(ctemp::Float64, response::PupaMortalityAbiodun, ::Float64)

Return `q_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Pupa. Source: Abiodun et al (2016).
"""
function get_temperature_response(ctemp::Float64, response::PupaMortalityAbiodun, ::Float64)
    _, qP = _get_coefficients(ctemp, response)
    return exp(-1 / qP)
end

"""
Data for `μ_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Pupa. Source: Abiodun et al (2016).
"""
mutable struct AdultMortalityAbiodun <: TemperatureResponse
    a::Float64
    b::Float64
    c::Float64
end

"""
    get_temperature_response(ctemp::Float64, response::AdultMortalityAbiodun, ::Float64)

Return `q_temperature_response`. `Species`: AnophelesGambiae, `LifeStage`: Male, Female. Source: Abiodun et al (2016).
"""
function get_temperature_response(
    ctemp::Float64,
    response::AdultMortalityAbiodun,
    ::Float64,
)
    # TODO: Update M & F to calculate lower mort rate/longer life for F
    return 1 / (-response.a + response.b * ctemp - response.c * ctemp^2)
end
