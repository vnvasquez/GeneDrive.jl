
# life stages/genes not specified retain default transition rate of zero
migration_data_five_stage_wolbachia_two_node = Dict(
    ("Male", "WW") => Dict((:Cuenca, :Guayaquil) => 0.002,
                           (:Guayaquil, :Cuenca) => 0.002),
    ("Male", "ww") => Dict((:Cuenca, :Guayaquil) => 0.002,
                            (:Guayaquil, :Cuenca) => 0.002),
    ("Female", "WW") => Dict((:Cuenca, :Guayaquil) => 0.002,
                            (:Guayaquil, :Cuenca) => 0.002),
    ("Female", "ww") => Dict((:Cuenca, :Guayaquil) => 0.002,
                            (:Guayaquil, :Cuenca) => 0.002),)
