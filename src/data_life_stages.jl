########################################
#          Density Dependence          #
########################################

dens_lin() = Density(LinearDensity, 1.0);
dens_log() = Density(LogisticDensity, 1.0);
dens_none() = Density(NoDensity, 1.0);

########################################
#   AedesAegypti Temperature Response  #
########################################

# Rossi et al 2014 (& Pelotti et al 2011 on pupa duration)
eggmortalityrossi = EggMortalityRossi(0.0731, 0.0595)
eggdurationrossi = EggDurationRossi(0.00764, 273.0, 40.55, 13094.10, 92.501, 28169.2)

larvamortalityrossi = LarvaMortalityRossi(0.0143, 0.00189)
larvadurationrossi = LarvaDurationRossi(0.00219, 273.0, 25.21, 7514.34)

pupamortalityrossi = PupaMortalityRossi(0.0143, 0.00189)
pupadurationrossi = PupaDurationRossi(0.027, -1.7, 27.7)

adultmortalityrossi = AdultMortalityRossi(0.053, 0.081, 23.0, 6.375 * 10^(-4))

function instantiate_stages_rossi()
    return DataStructures.OrderedDict(
        Egg => Stage{Egg}(eggmortalityrossi, eggdurationrossi, 2, dens_none(), nothing, 0),
        Larva =>
            Stage{Larva}(larvamortalityrossi, larvadurationrossi, 3, dens_log(), Egg, 0),
        Pupa =>
            Stage{Pupa}(pupamortalityrossi, pupadurationrossi, 2, dens_none(), Larva, 0),
        Male => Stage{Male}(adultmortalityrossi, 1, dens_none(), Pupa, 0),
        Female => Stage{Female}(adultmortalityrossi, 1, dens_none(), Pupa, 500),
    )
end

stages_rossi() = instantiate_stages_rossi();

# El Moustaid et al 2019
eggdurationmoustaid = EggDurationMoustaid(87.1722, 1.1575, 0.0136)
eggmortalitymoustaid = EggMortalityMoustaid(30.33, 9.395, 3.6576, 110.0014, 0.00077)

larvadurationmoustaid = LarvaDurationMoustaid(29.0889, 1.9571, 0.03357)
larvamortalitymoustaid = LarvaMortalityMoustaid(36.79, 2.354, 2.01938, 73.1928, 0.00046)

pupadurationmoustaid = PupaDurationMoustaid(20.5074, 1.0954, 0.0153)
pupamortalitymoustaid = PupaMortalityMoustaid(37.4679, 9.96278, 1.0, 38.0, 0.03487954)

adultmortalitymoustaid = AdultMortalityMoustaid(37.73, 9.16, -0.149)

function instantiate_stages_moustaid()
    return DataStructures.OrderedDict(
        Egg => Stage{Egg}(
            eggmortalitymoustaid,
            eggdurationmoustaid,
            2,
            dens_none(),
            nothing,
            0,
        ),
        Larva => Stage{Larva}(
            larvamortalitymoustaid,
            larvadurationmoustaid,
            3,
            dens_log(),
            Egg,
            0,
        ),
        Pupa => Stage{Pupa}(
            pupamortalitymoustaid,
            pupadurationmoustaid,
            2,
            dens_none(),
            Larva,
            0,
        ),
        Male => Stage{Male}(adultmortalitymoustaid, 1, dens_none(), Pupa, 0),
        Female => Stage{Female}(adultmortalitymoustaid, 1, dens_none(), Pupa, 500),
    )
end

stages_moustaid() = instantiate_stages_moustaid();

########################################
# AnophelesGambiae Temperature Response#
########################################

# Abiodun et al 2016
eggdurationabiodun = EggDurationAbiodun(6.0, 1.011, 20.212, 1.0, 12.096, 4.839, 1.0)
eggmortalityabiodun = EggMortalityAbiodun(6.0, 1.011, 20.212, 1.0, 12.096, 4.839, 1.0)

larvadurationabiodun =
    LarvaDurationAbiodun(6.0, 8.130, 13.794, 1.0, 12.096, 4.839, 1.0, 1.011, 20.212)
larvamortalityabiodun =
    LarvaMortalityAbiodun(6.0, 8.130, 13.794, 1.0, 12.096, 4.839, 1.0, 1.011, 20.212)

pupadurationabiodun = PupaDurationAbiodun(
    6.0,
    8.560,
    20.654,
    1.0,
    19.759,
    6.827,
    1.0,
    8.130,
    13.794,
    12.096,
    4.839,
)
pupamortalityabiodun = PupaMortalityAbiodun(
    6.0,
    8.560,
    20.654,
    1.0,
    19.759,
    6.827,
    1.0,
    8.130,
    13.794,
    12.096,
    4.839,
)

adultmortalityabiodun = AdultMortalityAbiodun(4.4, 1.31, 0.03)

function instantiate_stages_abiodun()
    return DataStructures.OrderedDict(
        Egg =>
            Stage{Egg}(eggmortalityabiodun, eggdurationabiodun, 2, dens_none(), nothing, 0),
        Larva => Stage{Larva}(
            larvamortalityabiodun,
            larvadurationabiodun,
            3,
            dens_log(),
            Egg,
            0,
        ),
        Pupa => Stage{Pupa}(
            pupamortalityabiodun,
            pupadurationabiodun,
            2,
            dens_none(),
            Larva,
            0,
        ),
        Male => Stage{Male}(adultmortalityabiodun, 1, dens_none(), Pupa, 0),
        Female => Stage{Female}(adultmortalityabiodun, 1, dens_none(), Pupa, 500),
    )
end

stages_abiodun() = instantiate_stages_abiodun();
