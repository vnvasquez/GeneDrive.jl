
################################################################################
#                              Genetic Structs                                 #
################################################################################

abstract type Genotype end

"""
        mutable struct Drive{G <: Genotype}
            genotype::Type{G}
            cube_slice::Array{Float64,2}
            s::Float64
            τ::Array{Float64,2}
            ϕ::Float64
            ξ_m::Float64
            ξ_f::Float64
            ω::Float64
            β::Float64
            η::Float64
            wildtype::Int64
            modified::Int64
        end

    Data for individual genotypes.

# Fields
- `genotype::Type{G}`: Single genotype.
- `cube_slice::Array{Float64,2}`: Offspring likelihoods for this genotype.
- `s::Float64`: Multiplicative fertility modifier.
- `τ::Array{Float64,2}`: Offspring viability.
- `ϕ::Float64`: Male to female emergence ratio.
- `ξ_m::Float64`: Male pupatory success.
- `ξ_f::Float64`: Female pupatory success.
- `ω::Float64`: Multiplicative adult mortality modifier.
- `β::Float64`: Female fecundity.
- `η::Float64`: Male mating fitness.
- `wildtype::Int64`: Boolean demarcates homozygous recessive.
- `modified::Int64`: Boolean demarcates homozygous modified.
"""
mutable struct Drive{G <: Genotype}
    genotype::Type{G}
    cube_slice::Array{Float64,2}
    #cube_slice::TemperatureResponse
    s::Float64
    τ::Array{Float64,2}
    #τ::TemperatureResponse
    ϕ::Float64
    ξ_m::Float64
    ξ_f::Float64
    ω::Float64
    β::Float64
    η::Float64
    wildtype::Int64
    modified::Int64
end

abstract type Construct end

"""
        mutable struct Genetics

            all_genotypes::Array{Drive{<:Genotype}}
            cube::Array{Float64, 3}
            S::Vector{Float64}
            Τ::Array{Float64,3}
            Φ::Vector{Float64}
            Ξ_m::Vector{Float64}
            Ξ_f::Vector{Float64}
            Ω::Vector{Float64}
            Β::Vector{Float64}
            Η::Vector{Float64}
            all_wildtypes::Vector{Int64}
            all_modified::Vector{Int64}

                function Genetics(all_genotypes::Array{Drive{<:Genotype}})

                    gN = length(all_genotypes)
                    cube = Array{Float64, 3}(undef, gN, gN, gN)
                    S = Vector{Float64}(undef, gN)
                    Τ = Array{Float64,3}(undef, gN, gN, gN)
                    Φ = Vector{Float64}(undef, gN)
                    Ξ_m = Vector{Float64}(undef, gN)
                    Ξ_f = Vector{Float64}(undef, gN)
                    Ω = Vector{Float64}(undef, gN)
                    Β = Vector{Float64}(undef, gN)
                    Η = Vector{Float64}(undef, gN)
                    all_wildtypes = Vector{Int64}(undef, gN)
                    all_modified = Vector{Int64}(undef, gN)

                    for (index, g) in enumerate(all_genotypes)
                        cube[:,:,index] = g.cube_slice
                        S[index] = g.s
                        Τ[:,:,index] = g.τ
                        Φ[index] = g.ϕ
                        Ξ_m[index] = g.ξ_m
                        Ξ_f[index] = g.ξ_f
                        Ω[index] = g.ω
                        Β[index] = g.β
                        Η[index] = g.η
                        all_wildtypes[index] = g.wildtype
                        all_modified[index] = g.modified
                    end

                    new(all_genotypes, cube, S, Τ, Φ, Ξ_m, Ξ_f, Ω, Β, Η,
                    all_wildtypes, all_modified)

                end

        end

        Data for all genotypes in a population.

# Fields
- `all_genotypes::Array{Drive{<:Genotype}}`: All genotypes in a population.
- `cube::Array{Float64, 3}`: Offspring likelihoods per genotype.
- `S::Vector{Float64}`: Multiplicative fertility modifier genedata_splitdriveper genotype, applied to oviposition.
- `Τ::Array{Float64,3}`: Offspring viability per genotype, applied to oviposition.
- `Φ::Vector{Float64}`: Male to female emergence ratio per genotype.
- `Ξ_m::Vector{Float64}`: Male pupatory success per genotype.
- `Ξ_f::Vector{Float64}`: Female pupatory success per genotype.
- `Ω::Vector{Float64}`: Multiplicative adult mortality modifier per genotype.
- `Β::Vector{Float64}`: Female fecundity per genotype (count of eggs laid).
- `Η::Vector{Float64}`: Male mating fitness per genotype.
- `all_wildtypes::Vector{Int64}`: Collect homozygous wildtype booleans.
- `all_modified::Vector{Int64}`: Collect homozygous modified booleans.
"""
mutable struct Genetics{C <: Construct}

    all_genotypes::Array{Drive{<:Genotype}}
    cube::Array{Float64, 3}
    #cube::Array{TemperatureResponse}
    S::Vector{Float64}
    Τ::Array{Float64,3}
    #Τ::Array{TemperatureResponse}
    Φ::Vector{Float64}
    Ξ_m::Vector{Float64}
    Ξ_f::Vector{Float64}
    Ω::Vector{Float64}
    Β::Vector{Float64}
    Η::Vector{Float64}
    all_wildtypes::Vector{Int64}
    all_modified::Vector{Int64}

        function Genetics(::Type{C}, all_genotypes::Array{Drive{<:Genotype}}) where {C <: Construct}

            gN = length(all_genotypes)
            cube = Array{Float64, 3}(undef, gN, gN, gN)
            S = Vector{Float64}(undef, gN)
            Τ = Array{Float64,3}(undef, gN, gN, gN)
            Φ = Vector{Float64}(undef, gN)
            Ξ_m = Vector{Float64}(undef, gN)
            Ξ_f = Vector{Float64}(undef, gN)
            Ω = Vector{Float64}(undef, gN)
            Β = Vector{Float64}(undef, gN)
            Η = Vector{Float64}(undef, gN)
            all_wildtypes = Vector{Int64}(undef, gN)
            all_modified = Vector{Int64}(undef, gN)

            for (index, g) in enumerate(all_genotypes)
                cube[:,:,index] = g.cube_slice
                S[index] = g.s
                Τ[:,:,index] = g.τ
                Φ[index] = g.ϕ
                Ξ_m[index] = g.ξ_m
                Ξ_f[index] = g.ξ_f
                Ω[index] = g.ω
                Β[index] = g.β
                Η[index] = g.η
                all_wildtypes[index] = g.wildtype
                all_modified[index] = g.modified
            end

            new{C}(all_genotypes, cube, S, Τ, Φ, Ξ_m, Ξ_f, Ω, Β, Η,
            all_wildtypes, all_modified)

        end

end


########################################
#               Wolbachia              #
########################################
struct Wolbachia <: Construct end
struct WW <: Genotype end
struct ww <: Genotype end # wildtype

########################################
#                  RIDL                #
########################################
struct RIDL <: Construct end
#struct WW <: Genotype end # wildtype
#struct WR <: Genotype end
#struct RR <: Genotype end # homozygous modified

########################################
#                 MCR                  #
########################################
struct MCR <: Construct end
struct HH <: Genotype end # homozygous modified
struct Hh <: Genotype end # heterozygous
struct HR <: Genotype end
struct hh <: Genotype end # wildtype
struct hR <: Genotype end
struct RR <: Genotype end

########################################
#         Single Locus Homing          #
########################################
struct SLHoming <: Construct end
# struct WW <: Genotype end # wildtype
struct WH <: Genotype end
struct WR <: Genotype end
struct WB <: Genotype end
# struct HH <: Genotype end # homozygous modified
# struct HR <: Genotype end
struct HB <: Genotype end
# struct RR <: Genotype end
struct RB <: Genotype end
struct BB <: Genotype end

########################################
#               Split Drive            #
########################################
struct SplitDrive <: Construct end
struct WWWW <: Genotype end # wildtype
struct WWWH <: Genotype end
struct WWWR <: Genotype end
struct WWWB <: Genotype end
struct WWHH <: Genotype end
struct WWHR <: Genotype end
struct WWHB <: Genotype end
struct WWRR <: Genotype end
struct WWRB <: Genotype end
struct WWBB <: Genotype end
struct WCWW <: Genotype end
struct WCWH <: Genotype end
struct WCWR <: Genotype end
struct WCWB <: Genotype end
struct WCHH <: Genotype end
struct WCHR <: Genotype end
struct WCHB <: Genotype end
struct WCRR <: Genotype end
struct WCRB <: Genotype end
struct WCBB <: Genotype end
struct CCWW <: Genotype end
struct CCWH <: Genotype end
struct CCWR <: Genotype end
struct CCWB <: Genotype end
struct CCHH <: Genotype end
struct CCHR <: Genotype end
struct CCHB <: Genotype end
struct CCRR <: Genotype end
struct CCRB <: Genotype end
struct CCBB <: Genotype end
