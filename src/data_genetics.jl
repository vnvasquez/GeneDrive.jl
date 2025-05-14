# Data source: https://cran.r-project.org/web/packages/MGDrivE/MGDrivE.pdf
########################################
#              Mendelian               #
########################################

mend1 = [1.0 0.50 0.0; 0.5 0.25 0.0; 0.0 0.0 0.0]
mend2 = [0.0 0.5 1.0; 0.5 0.5 0.5; 1.0 0.5 0.0]
mend3 = [0.0 0.00 0.0; 0.0 0.25 0.5; 0.0 0.50 1.0]

drives_mendelian = [
    Drive(AA, mend1, 1.0, ones(3, 3), 0.5, 1.0, 1.0, 0.0, 63.0, 1.0, 0, 1),
    Drive(Aa, mend2, 1.0, ones(3, 3), 0.5, 1.0, 1.0, 0.0, 63.0, 1.0, 0, 0),
    Drive(aa, mend3, 1.0, ones(3, 3), 0.5, 1.0, 1.0, 0.0, 63.0, 1.0, 1, 0),
]

"""
    genetics_mendelian()

Return Mendelian inheritance. 
"""
genetics_mendelian() = Genetics(Mendelian, drives_mendelian)

########################################
#                  RIDL                #
########################################

ridl1 = [1.0 0.5 0.0; 0.5 0.25 0.0; 0.0 0.0 0.0]
ridl2 = [0.0 0.5 1.0; 0.5 0.5 0.5; 1.0 0.5 0.0]
ridl3 = [0.0 0.0 0.0; 0.0 0.25 0.5; 0.0 0.50 1.0]

drives_ridl = [ # NB: Midrange estimates for ω and η
    Drive(WW, ridl1, 1.0, ones(3, 3), 0.5, 1.0, 1.0, 0.0, 63.0, 1.0, 1, 0),
    Drive(WR, ridl2, 1.0, ones(3, 3), 0.5, 0.0, 0.0, 0.18, 63.0, 0.031, 0, 0),
    Drive(RR, ridl3, 1.0, ones(3, 3), 0.5, 0.0, 0.0, 0.18, 63.0, 0.031, 0, 1),
]

"""
    genetics_ridl()

Return Release of Insects carrying a Dominant Lethal gene (RIDL) inheritance. 
"""
genetics_ridl() = Genetics(RIDL, drives_ridl)

########################################
#               Wolbachia              #
########################################

wolb1 = [1.0 1.0; 1.0 0.0]
wolb2 = [0.0 0.0; 0.0 1.0]

tau1_wolb = [1.0 1.0; 0.0 1.0]
tau2_wolb = [1.0 1.0; 0.0 1.0]

drives_wolbachia = [ # NB: Midrange estimate for ω
    Drive(WW, wolb1, 1.0, tau1_wolb, 0.5, 1.0, 1.0, 0.1, 63.0, 1.0, 0, 1),
    Drive(ww, wolb2, 1.0, tau2_wolb, 0.5, 1.0, 1.0, 0.0, 63.0, 1.0, 1, 0),
]

"""
    genetics_wolbachia()

Return Wolbachia inheritance. 
"""
genetics_wolbachia() = Genetics(Wolbachia, drives_wolbachia)

########################################
#              HGD/MCR                 #
########################################

mcr1 = [
    1.0 1.0 0.50 0 0 0
    1.0 1.0 0.50 0 0 0
    0.5 0.5 0.25 0 0 0
    0.0 0.0 0.00 0 0 0
    0.0 0.0 0.00 0 0 0
    0.0 0.0 0.00 0 0 0
]
mcr2 = [
    0.0 0.0 0.00 1.0 0.50 0
    0.0 0.0 0.00 1.0 0.50 0
    0.0 0.0 0.00 0.5 0.25 0
    1.0 1.0 0.50 0.0 0.00 0
    0.5 0.5 0.25 0.0 0.00 0
    0.0 0.0 0.00 0.0 0.00 0
]
mcr3 = [
    0.0 0.0 0.50 0 0.50 1.0
    0.0 0.0 0.50 0 0.50 1.0
    0.5 0.5 0.50 0 0.25 0.5
    0.0 0.0 0.00 0 0.00 0.0
    0.5 0.5 0.25 0 0.00 0.0
    1.0 1.0 0.50 0 0.00 0.0
]
mcr4 = [
    0 0 0 0.0 0.00 0
    0 0 0 0.0 0.00 0
    0 0 0 0.0 0.00 0
    0 0 0 1.0 0.50 0
    0 0 0 0.5 0.25 0
    0 0 0 0.0 0.00 0
]
mcr5 = [
    0 0 0 0.0 0.00 0
    0 0 0 0.0 0.00 0
    0 0 0.00 0.5 0.25 0.0
    0 0 0.50 0.0 0.50 1.0
    0 0 0.25 0.5 0.50 0.5
    0 0 0.00 1.0 0.50 0.0
]
mcr6 = [
    0 0 0 0.0 0.00 0
    0 0 0 0.0 0.00 0
    0 0 0.25 0 0.25 0.5
    0 0 0.00 0 0.00 0.0
    0 0 0.25 0 0.25 0.5
    0 0 0.50 0 0.50 1.0
]

drives_mcr = [ # NB: ∉ costs
    Drive(HH, mcr1, 1.0, ones(6, 6), 0.5, 1.0, 1.0, 0.0, 63.0, 1.0, 0, 1),
    Drive(Hh, mcr2, 1.0, ones(6, 6), 0.5, 1.0, 1.0, 0.0, 63.0, 1.0, 0, 0),
    Drive(HR, mcr3, 1.0, ones(6, 6), 0.5, 1.0, 1.0, 0.0, 63.0, 1.0, 0, 0),
    Drive(hh, mcr4, 1.0, ones(6, 6), 0.5, 1.0, 1.0, 0.0, 63.0, 1.0, 1, 0),
    Drive(hR, mcr5, 1.0, ones(6, 6), 0.5, 1.0, 1.0, 0.0, 63.0, 1.0, 0, 0),
    Drive(RR, mcr6, 1.0, ones(6, 6), 0.5, 1.0, 1.0, 0.0, 63.0, 1.0, 0, 0),
]


"""
    genetics_mcr()

Return Mutagenic Chain Reaction (MCR) inheritance. 
"""
genetics_mcr() = Genetics(MCR, drives_mcr)
