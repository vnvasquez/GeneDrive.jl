# TODO: consider exporting abstract types later
abstract type TempResponse end

function get_temperature_response(ctemp::Float64, response::TempResponse)
    return get_temperature_response(ctemp::Float64, response::TempResponse, 0.0)
end

"""
    Data for model without temperature/temperature response. Applies to all species,
    all life stages.
"""
mutable struct NoResponse <: TempResponse end

"""
    Function for model without temperature/temperature response. Applies to all species,
    all life stages.
"""
get_temperature_response(::Float64, ::NoResponse, ::Float64) = 0.0

"""
    mutable struct EggMortalityRossi <: TempResponse
        a::Float64
        b::Float64
    end

    Data for temperature sensitive mortality. Applies to AedesAegypti, egg stage.
    Source: Rossi et al (2014).

# Arguments
- `a::Float64`: Death rate at low temperatures, validation range: `(0, nothing)`
- `b::Float64`: Influence factor of temperature, validation range: `(0, nothing)`
"""
mutable struct EggMortalityRossi <: TempResponse
    a::Float64
    b::Float64
end

"""
    Function for temperature sensitive mortality. Applies to AedesAegypti,
    egg stage. Source: Rossi et al (2014).
"""
function get_temperature_response(ctemp::Float64, response::EggMortalityRossi)
    return response.a * exp(response.b * ctemp)
end

"""
    Data for temperature sensitive duration. Applies to AedesAegypti, egg stage.
    Source: Rossi et al (2014).
"""
mutable struct EggDurationRossi <: TempResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
end

"""
    Function for temperature sensitive duration. Applies to AedesAegypti, egg stage.
    Source: Rossi et al (2014).
"""
function get_temperature_response(ctemp::Float64, response::EggDurationRossi, ::Float64)
    return 1 / (response.a*(ctemp+response.b) *
        (exp(response.c - (response.d/(ctemp+response.b))) /
        1+exp(response.e - (response.f /(ctemp + response.b)))))
end

"""
    Data for temperature sensitive mortality. Applies to AedesAegypti, larval stage.
    Source: Rossi et al (2014).
"""
mutable struct LarvaMortalityRossi <: TempResponse
    a::Float64
    b::Float64
    # TODO: why temp dependent deathrate unused in source
    # c:: Float64 # π0 = 4.03*10^(-6)
end

"""
    Function for temperature sensitive mortality. Applies to AedesAegypti, larval stage.
    Source: Rossi et al (2014).
"""
function get_temperature_response(ctemp::Float64, response::LarvaMortalityRossi, ::Float64)
    return response.a * exp(response.b * ctemp^2)
end

"""
    Data for temperature sensitive duration. Applies to AedesAegypti, larval stage.
    Source: Rossi et al (2014).
"""
mutable struct LarvaDurationRossi <: TempResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
end

"""
    Function for temperature sensitive duration. Applies to AedesAegypti, larval stage.
    Source: Rossi et al (2014).
"""
function get_temperature_response(ctemp::Float64, response::LarvaDurationRossi, ::Float64)
    return 1 / (response.a *(ctemp + response.b) *
        exp(response.c - (response.d / (ctemp + response.b))))
end


"""
    Data for temperature sensitive mortality. Applies to AedesAegypti, pupal stage.
    Source: Rossi et al (2014).
"""
mutable struct PupaMortalityRossi <: TempResponse
    a::Float64
    b::Float64
    # TODO: why temp dependent deathrate unused in source
    # c:: Float64 # π0 = 4.03*10^(-6)
end

"""
    Function for temperature sensitive mortality. Applies to AedesAegypti, pupal stage.
    Source: Rossi et al (2014).
"""
function get_temperature_response(ctemp::Float64, response::PupaMortalityRossi, ::Float64)
    return response.a * exp(response.b * ctemp^2)
end

"""
    Data for temperature sensitive duration. Applies to AedesAegypti, pupal stage.
    Source: Rossi et al (2014).
"""
mutable struct PupaDurationRossi <: TempResponse
    a::Float64
    b::Float64
    c::Float64
end

"""
    Function for temperature sensitive duration. Applies to AedesAegypti, pupal stage.
    Source: Rossi et al (2014).
"""
function get_temperature_response(ctemp::Float64, response::PupaDurationRossi, ::Float64)
    return response.a * ctemp^2 + response.b * ctemp + response.c
end

"""
    Data for temperature sensitive mortality. Applies to AedesAegypti, adult stage.
    Source: Rossi et al (2014).
"""
mutable struct AdultMortalityRossi <: TempResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
end

"""
    Function for temperature sensitive mortality. Applies to AedesAegypti, adult stage.
    Source: Rossi et al (2014).
"""
function get_temperature_response(ctemp::Float64, response::AdultMortalityRossi, ::Float64)
    return response.a * exp( -response.b *
        (ctemp - response.c)^2) + response.d * (ctemp - response.c)^2
end

"""
    Data for temperature sensitive duration. Applies to AedesAegypti, egg stage.
    Source: El Moustaid et al (2019).
"""
mutable struct EggDurationMoustaid <: TempResponse
    a::Float64
    b::Float64
    c::Float64
end

"""
    Function for temperature sensitive duration. Applies to AedesAegypti, egg stage.
    Source: El Moustaid et al (2019).
"""
function get_temperature_response(ctemp::Float64, response::EggDurationMoustaid, ::Float64)
    # TODO: why q = 1/duration in original implementation (applies to all juv), and is ~65 days reasonable?
    return 1/(response.a - (response.b * ctemp) + response.c * ctemp^2)
end

"""
    Data for temperature sensitive mortality. Applies to AedesAegypti, egg stage.
    Source: El Moustaid et al (2019).
"""
mutable struct EggMortalityMoustaid <: TempResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
end

"""
    Function for temperature sensitive mortality. Applies to AedesAegypti, egg stage.
    Source: El Moustaid et al (2019).
"""
function get_temperature_response(ctemp::Float64, response::EggMortalityMoustaid, duration::Float64)
    survival = 0.0
    if ctemp > response.a || ctemp < response.b
        survival = 0.0
    elseif ctemp > response.a && ctemp <= response.a
        survival = -ctemp * response.c + response.d
    else
        survival = response.e * ctemp * (ctemp - response.b) * ((response.a - ctemp)^0.5)
    end
    return (1 - survival) / (1/duration)
end

"""
    Data for temperature sensitive duration. Applies to AedesAegypti, larval stage.
    Source: El Moustaid et al (2019).
"""
mutable struct LarvaDurationMoustaid <: TempResponse
    a::Float64
    b::Float64
    c::Float64
end

"""
    Function for temperature sensitive duration. Applies to AedesAegypti, larval stage.
    Source: El Moustaid et al (2019).
"""
function get_temperature_response(ctemp::Float64, response::LarvaDurationMoustaid, ::Float64)
    return 1 / (response.a - (response.b * ctemp) + response.c * ctemp^2)
end

"""
    Data for temperature sensitive mortality. Applies to AedesAegypti, larval stage.
    Source: El Moustaid et al (2019).
"""
mutable struct LarvaMortalityMoustaid <: TempResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
end

"""
    Function for temperature sensitive mortality. Applies to AedesAegypti, larval stage.
    Source: El Moustaid et al (2019).
"""
function get_temperature_response(ctemp::Float64, response::LarvaMortalityMoustaid, duration::Float64)
    survival = 0.0
    if ctemp > response.a || ctemp < response.b
        survival = 0.0
    elseif ctemp > response.a && ctemp <= response.a
        survival = -ctemp * response.c + response.d
    else
        survival = response.e * ctemp * (ctemp - response.b) * ((response.a - ctemp)^0.5)
    end
    return (1 - survival) / (1/duration)
end

"""
    Data for temperature sensitive duration. Applies to AedesAegypti, pupal stage.
    Source: El Moustaid et al (2019).
"""
mutable struct PupaDurationMoustaid <: TempResponse
    a::Float64
    b::Float64
    c::Float64
end

"""
    Function for temperature sensitive duration. Applies to AedesAegypti, pupal stage.
    Source: El Moustaid et al (2019).
"""
function get_temperature_response(ctemp::Float64, response::PupaDurationMoustaid, ::Float64)
    return  1 / (response.a - (response.b * ctemp) + response.c * ctemp^2)
end

"""
    Data for temperature sensitive mortality. Applies to AedesAegypti, pupal stage.
    Source: El Moustaid et al (2019).
"""
mutable struct PupaMortalityMoustaid <: TempResponse
    a::Float64
    b::Float64
    c::Float64 # El Moustaid personal comm, lit value incorrect
    d::Float64
    e::Float64
end

"""
    Function for temperature sensitive mortality. Applies to AedesAegypti, pupal stage.
    Source: El Moustaid et al (2019).
"""
function get_temperature_response(ctemp::Float64, response::PupaMortalityMoustaid, duration::Float64)
    survival = 0.0
    if ctemp > response.a || ctemp < response.b
        survival = 0.0
    elseif ctemp > response.a && ctemp <= response.a
        survival = -ctemp * response.c + response.d
    else
        survival = min(1,response.e * ctemp * (ctemp - response.b) * ((response.a - ctemp)^0.5))
    end
    return (1 - survival) / (1 / duration)
end

"""
    Data for temperature sensitive mortality. Applies to AedesAegypti, adult stage.
    Source: El Moustaid et al (2019).
"""
mutable struct AdultMortalityMoustaid <: TempResponse
    a::Float64 # El Moustaid personal comm, lit value incorrect
    b::Float64 # El Moustaid personal comm, lit value incorrect
    c::Float64 # El Moustaid personal comm: use "lf" param from Mordecai 2017, Table B
end


"""
    Function for temperature sensitive mortality. Applies to AedesAegypti, adult stage.
    Source: El Moustaid et al (2019).
"""
function get_temperature_response(ctemp::Float64, response::AdultMortalityMoustaid, ::Float64)

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
        μ = 1/survival
    end

    return μ
end


"""
    Data for temperature sensitive duration. Applies to AnophelesGambiae, egg stage.
    Source: Abiodun et al (2016).
"""
mutable struct EggDurationAbiodun <: TempResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
end

"""
    Data for temperature sensitive mortality. Applies to AnophelesGambiae, egg stage.
    Source: Abiodun et al (2016).
"""
mutable struct EggMortalityAbiodun <: TempResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
end

# Mortality and duration use the same coefficients
function _get_coefficients(ctemp::Float64, response::Union{EggDurationAbiodun, EggMortalityAbiodun})
    ctemp += response.a
    q = response.b + response.c * (response.d + (ctemp/response.e)^response.f)^(-response.g)
    return q
end

"""
    Function for temperature sensitive duration. Applies to AnophelesGambiae, egg stage.
    Source: Abiodun et al (2016).
"""
function get_temperature_response(ctemp::Float64, response::EggDurationAbiodun, ::Float64)
    return _get_coefficients(ctemp, response)
end

"""
    Function for temperature sensitive duration. Applies to AnophelesGambiae, egg stage.
    Source: Abiodun et al (2016).
"""
function get_temperature_response(ctemp::Float64, response::EggMortalityAbiodun, ::Float64)
    q = _get_coefficients(ctemp, response)
    return exp(-1/q)
end

"""
    Data for temperature sensitive duration. Applies to AnophelesGambiae, larva stage.
    Source: Abiodun et al (2016).
"""
mutable struct LarvaDurationAbiodun <: TempResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
end

"""
    Data for temperature sensitive mortality. Applies to AnophelesGambiae, larva stage.
    Source: Abiodun et al (2016).
"""
mutable struct LarvaMortalityAbiodun <: TempResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
end

# Mortality and duration use the same coefficients
function _get_coefficients(ctemp::Float64, response::Union{LarvaDurationAbiodun, LarvaMortalityAbiodun})
    ctemp += response.a
    q = response.b + response.c * (response.d + (ctemp/response.e)^response.f)^(-response.g)
    return q
end

"""
    Function for temperature sensitive duration. Applies to AnophelesGambiae, larva stage.
    Source: Abiodun et al (2016).
"""
function get_temperature_response(ctemp::Float64, response::LarvaDurationAbiodun, egg_duration::Float64)
    q = _get_coefficients(ctemp, response)
    return q - egg_duration
end

"""
    Function for temperature sensitive duration. Applies to AnophelesGambiae, larva stage.
    Source: Abiodun et al (2016).
"""
function get_temperature_response(ctemp::Float64, response::LarvaMortalityAbiodun, ::Float64)
    q = _get_coefficients(ctemp, response)
    return exp(-1/q)
end

"""
    Data for temperature sensitive duration. Applies to AnophelesGambiae, pupa stage.
    Source: Abiodun et al (2016).
"""
mutable struct PupaDurationAbiodun <: TempResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
end

"""
    Data for temperature sensitive mortality. Applies to AnophelesGambiae, pupa stage.
    Source: Abiodun et al (2016).
"""
mutable struct PupaMortalityAbiodun <: TempResponse
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
end

# Mortality and duration use the same coefficients
function _get_coefficients(ctemp::Float64, response::Union{PupaDurationAbiodun, PupaMortalityAbiodun})
    ctemp += response.a
    q = response.b + response.c * (response.d + (ctemp/response.e)^response.f)^(-response.g)
    return q
end

"""
    Function for temperature sensitive duration. Applies to AnophelesGambiae, pupa stage.
    Source: Abiodun et al (2016).
"""
function get_temperature_response(ctemp::Float64, response::PupaDurationAbiodun, larva_duration::Float64)
    q = _get_coefficients(ctemp, response)
    return q - larva_duration
end

"""
    Function for temperature sensitive duration. Applies to AnophelesGambiae, pupa stage.
    Source: Abiodun et al (2016).
"""
function get_temperature_response(ctemp::Float64, response::PupaMortalityAbiodun, ::Float64)
    q = _get_coefficients(ctemp, response)
    return exp(-1/q)
end

"""
    Data for temperature sensitive mortality. Applies to AnophelesGambiae, pupa stage.
    Source: Abiodun et al (2016).
"""
mutable struct AdultMortalityAbiodun <: TempResponse
    a::Float64
    b::Float64
    c::Float64
end

"""
    Function for temperature sensitive duration. Applies to AnophelesGambiae, adult stage.
    Source: Abiodun et al (2016).
"""
function get_temperature_response(ctemp::Float64, response::AdultMortalityAbiodun, ::Float64)
    # TODO: Update M & F to calculate lower mort rate/longer life for F
    return 1 / (-response.a + response.b * ctemp - response.c * ctemp^2)
end


#=
function _get_coefficients(ctemp::Float64, input::Union{EggMortalityAbiodun, EggDurationAbiodun})
    # Tw = Ta + 2celsius (source: Table 3)
    ctemp += 6

    # mortality relies on duration; calc this first
    q_calc2 = 8.130 + 13.794 * (1 + (ctemp/12.096)^4.839)^(-1)

    # need duration from previous without altering arguments
    q_calc3 = 8.560 + 20.654 * (1 + (ctemp/19.759)^6.827)^(-1)

    return q_calc2, q_calc3
end

function get_temperature_response(ctemp::Float64, response::FakeResponseMortality, ::Float64)
    _, q_calc3 = _get_coefficients(ctemp, response)
    return exp(-1/q_calc3)
end

function get_temperature_response(ctemp::Float64, response::FakeResponseDuration, ::Float64)
    q_calc2, q_calc3 = _get_coefficients(ctemp, response)
    return q_calc3 - q_calc2
end
=#
