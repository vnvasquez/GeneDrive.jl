abstract type AbstractTempResponse end


get_temperature_response(ctemp::Nothing, ::AbstractTempResponse) = 0.0

"""
This is the data for the model from ref
"""
mutable struct ExampleResponse <: AbstractTempResponse
    a::Float64
    b::Float64
end

function get_temperature_response(ctemp::Float64, response::ExampleResponse)
    return ctemp*response.a + ctemp^2*response.b
end

