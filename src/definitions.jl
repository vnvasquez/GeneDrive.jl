################################################################################
#                                   Mappings                                   #
################################################################################

const life_stage_key_map =
    Dict("Egg" => Egg, "Larva" => Larva, "Pupa" => Pupa, "Male" => Male, "Female" => Female)

const genetics_key_map = Dict(
    "WW" => WW,
    "ww" => ww,
    # end Wolbachia
    "HH" => HH,
    "Hr" => Hh,
    "HR" => HR,
    "hh" => hh,
    "hR" => hR,
    "RR" => RR,
    ## end single-locus CRISPR homing ("MCR")
    "WR" => WR,
    # end RIDL (total = WW, WR, RR)
)
