################################################################################
#                               Organisms by Species                           #
################################################################################

abstract type Species end

"""
    mutable struct Organism{S <: Species}
        gene_data::Genetics
        all_stages::DataStructures.OrderedDict{Type{<:LifeStage}, Stage}
    end

Generic data container for species-specific genetic and life stage information.

# Arguments

  - `gene_data::Genetics`: Species and modification-specific genetic data.
  - `all_stages::DataStructures.OrderedDict{Type{<:LifeStage}, Stage}`: Dictionary of wildtype species-specific life stage data.
"""
mutable struct Organism{S <: Species}
    gene_data::Genetics
    all_stages::DataStructures.OrderedDict{Type{<:LifeStage}, Stage}
end

"""
    struct AedesAegypti <: Species end

Data for Aedes Aegypti mosquito. Disease vector.
"""
struct AedesAegypti <: Species end

"""
    struct AedesAlbopictus <: Species end

Data for Aedes Albopictus mosquito. Disease vector.
"""
struct AedesAlbopictus <: Species end

"""
    struct AnophelesGambiae <: Species end

Data for Anopheles Gambiae mosquito. Disease vector.
"""
struct AnophelesGambiae <: Species end

"""
    struct AnophelesArabiensis <: Species end

Data for Anopheles Arabiensis mosquito. Disease vector.
"""
struct AnophelesArabiensis <: Species end
