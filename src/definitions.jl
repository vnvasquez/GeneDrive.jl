################################################################################
#                                   Mappings                                   #
################################################################################



const life_stage_key_map = 

        Dict(
            "Egg" => Egg, 
            "Larva" => Larva, 
            "Pupa" => Pupa, 
            "Male" => Male, 
            "Female" => Female)

const genetics_key_map = 
    Dict(
        "WW" => WW,
        "ww" => ww,
        # end Wolbachia
        "HH" => HH,
        "Hh" => Hh,
        "HR" => HR,
        "hh" => hh,
        "hR" => hR,
        "RR" => RR,
        # end single-locus homing gene drive (HGD/MCR)
        "WR" => WR,
        # end RIDL (total = WW, WR, RR)
        "AA" => AA,
        "Aa" => Aa,
        "aa" => aa,
        # end Mendelian
        )

const MALE_FEMALE_RELEASE_FRACTION = 0.5






