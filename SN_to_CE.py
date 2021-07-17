# supernova.py

import os
import numpy as np
import mesa_reader as mr
import matplotlib.pyplot as plt

path = os.getcwd()

print('ATTENTION: have you input the values for SN kick from ComBinE?')

# constants
G = 6.67428e-8  # in g^-1 cm^3 s^-2
msol = 1.9892e33  # solar mass (g)
rsol = 6.9598e10  # solar radius (cm)
c = 29979245800  # speed of light in cm/s
G_sol = 3.92839e8  # in rsol**3 msol**-1 yr**-2

hist1 = mr.MesaData('../LOGS1/history.data')
histbi = mr.MesaData('../binary_history.data')
hist2 = mr.MesaData('../LOGS2/history.data')


model = hist1.model_number[-1]  # last model of history.data should be SN, use as index for binary history
model = model - 1

age = hist1.star_age[model]  # time of SN

M1 = msol * histbi.star_1_mass[model]  # pre-supernova mass of primary (star_1_mass)
M2 = msol * histbi.star_2_mass[model]  # mass of secondary (star_2_mass)
print('M1 = ', histbi.star_1_mass[model], 'Msun')
print('M2 = ', histbi.star_2_mass[model], 'Msun')
# M1 = 64.0284 * msol  # ComBinE
# M2 = 38.1828 * msol  # ComBinE
M0 = M1 + M2  # total pre-supernova mass

# get remnant mass
McoreCO = msol * hist1.c_core_mass[model]  # c_core_mass
print('mass of CO core in Msun: ', McoreCO/msol)
MlayerHe = msol * hist1.he_rich_layer_mass[model]  # he_rich_layer
print('mass of helium rich layer in Msun: ', MlayerHe/msol)

M1rem = 0.8 * (McoreCO + 0.8 * MlayerHe)  # ComBinE recipe from line 1239 in evolution.cpp, see Kruckow 2018 2.2.6
print(M1rem)
# M1rem = 23.3292 * msol  # ComBinE
M1remsol = M1rem / msol

DeltaM = M1 - M1rem  # mass ejected from SN
a0 = rsol * histbi.binary_separation[model]  # pre-SN binary separation (binary_separation)
# a0 = 3176.67 * rsol  # ComBinE
#e0 = histbi.eccentricity[model]  # pre-SN eccentricity (eccentricity) = 0
e0 = 0.0
print('initial eccentricity:', e0)
r = a0 * (1 + (e0 ** 2) / 2)
print('kick velocity in km/s: ')
w_in = float(input())  # km/s
w = w_in * (10 ** 5)  # kick velocity (from ComBinE screen output) in cm/s
theta = 90  # kick angle (from ComBinE screen output) in degrees
# np.cos needs input in rad
theta_rad = theta * np.pi / 180.0
phi = 90  # from ComBinE output
phi_rad = phi*np.pi/180.0
print('theta_rad = ', theta_rad)
print('cos(theta_rad) = ', np.cos(theta_rad))

print('DeltaM/M0: ', DeltaM / M0)



# test whether system becomes unbound by SN, according to J.G. Hills (1983)
# if more than half of system mass is lost by SN
if DeltaM / M0 >= (1 + (e0 ** 2)) / 2:
    print('System is destroyed by the SN')

else:
    vrel = (G * M0 / r) ** (1 / 2)

    print('vrel: ', vrel)
    print('w/vrel: ', w / vrel)
    print('a0/r: ', a0 / r)

    massfrac = DeltaM / M0

    # use 16.33 from Tauris & van den Heuvel 2006
    if theta == 90:
        afrac = ((1 - massfrac) / (1 - 2 * massfrac - ((w / vrel) ** 2)))  # a/a0
    else:
        afrac = ((1 - massfrac) / (1 - 2 * massfrac - ((w / vrel) ** 2) - 2 * np.cos(theta_rad) * (w / vrel)))  # a/a0


    print('a/a0 = ', afrac)

    a = afrac * a0
    a_rsol = a / rsol

    # calculate post-SN eccentricity with eq. 16.34 and 16.35
    mu = M1rem * M2 / (M1rem + M2)  # reduced mass of system

    if theta == 90:
        Lorb = r * mu * np.sqrt((vrel + (w * np.sin(phi_rad)) ** 2))
    else:
        Lorb = r * mu * np.sqrt((vrel + w * np.cos(theta_rad)) ** 2 + (w * np.sin(theta_rad) * np.sin(phi_rad)) ** 2)

    print('Lorb = ', Lorb)

    Eorb = - (G * M1rem * M2 / r) + (mu * (vrel ** 2) / 2)
    print('Eorb = ', Eorb)

    e_final = np.sqrt(1 + (2 * Eorb * (Lorb ** 2)) / (mu * (G ** 2) * (M1rem ** 2) * (M2 ** 2)))

    # re-circularize orbit immediately (not realistic, but result should not be significantly different)
    if w < 0.1:  # no kick: no change in orbital separation
        recircul = a0
        recircul_rsol = a0/rsol
    else:
        recircul = (1 - e_final**2)*a
        recircul_rsol = (1 - e_final**2)*a_rsol



    print('SN result for: ', path)
    print('system age: ', age)
    print('remnant mass in Msun: ', M1remsol)
    print('pre-SN masses in Msun: ', M1 / msol, ' ', M2 / msol)
    print('pre-Sn binary separation in Rsun: ', a0 / rsol)
    P_orb = 2 * np.pi * np.sqrt((a0 ** 3) / (G * (M1 + M2)))  # in seconds
    print('orbital period in days pre-SN: ', P_orb / 86400)
    #print('System survives SN with binary separation in Rsun: ', a_rsol)
    print('remnant mass in g = ', M1rem)
    # calculate Schwarzschild radius
    r_schw = 2 * G * M1rem / (c**2)
    print('corresponding Schwarzschild radius in cm: ', r_schw)
    print('Eccentricity of system: ', e_final)
    print('binary separation in Rsun before orbit is re-circularized: ', a_rsol)
    P_orb = 2 * np.pi * np.sqrt((a ** 3) / (G * (M1rem + M2)))  # in seconds
    print('orbital period in days before recircularisation: ', P_orb / 86400)
    print('binary separation in Rsun after orbit is re-circularized: ', recircul_rsol)
    P_orb = 2 * np.pi * np.sqrt((recircul ** 3) / (G * (M1rem + M2)))  # in seconds
    print('orbital period in days after recircularisation: ', P_orb / 86400)

    # radiustest.py


    index = 0  # indexing of mesa reader lists
    lastmodel = hist2.model_number[-1]  # last model number, where for loop ends
    end = lastmodel - 1  # last index (starts with 0)
    print('last model of secondary:', lastmodel)

    # import radius and mass lists
    R2 = hist2.radius  # Rsun
    M2 = hist2.star_mass
    age = hist2.star_age
    a0 = recircul_rsol  # compare with Rsun value

    # make list with Roche lobe radii
    R_L = np.zeros(lastmodel)

    for i in range(index, end): # start at 0 for Roche lobe calculation
        # Roche lobe equation
        q = M2[i] / M1remsol
        # calculate Roche lobe radius in Rsun from Eggleton 1983
        RL = 0.49 * q ** (2 / 3) * a0 / (0.6 * q ** (2 / 3) + np.log(1 + q ** (1 / 3)))
        R_L[i] = RL
        if R2[i] >= RL:
            print('secondary star expands to Roche lobe radius')
            print('model number where this happens (for use in common_envelope.py): ', hist2.model_number[i])
            print('time from ZAMS when CE happens: ', hist2.age[i])
            print('radius of secondary: ', R2[i])
            print('mass of secondary in Msun: ', M2[i])
            RL_last = RL
            crossover = i
            break

    # filling Roche lobe radius array with last calculated value for better indication
    # mention below figure
    for j in range(crossover, lastmodel):
        R_L[j] = RL_last

    plt.xlabel('t (Myr)')
    plt.ylabel('$R / \mathrm{R_{\odot}}$')
    plt.tick_params(axis='both', direction='in', which='both', top=True, right=True)

    # mass ratio
    ratio = M2[crossover] / M1remsol
    print('mass ratio at RLO: ', ratio)

    # second largest convective zone
    mix2_top = hist2.conv_mx2_top
    mix2_bot = hist2.conv_mx2_bot
    # second largest convective zone mass delta
    delta_mix = mix2_top[crossover] - mix2_bot[crossover]

    # there may be more than one convective region in the envelope
    # this is just a simplified test
    print('CE initiated if convective envelope > 0.1 in mass coordinate, largest convective region in envelope has this value in m/mstar: ', delta_mix)


    # common_envelope.py

    model = crossover  # result from radiustest.py
    bind_g = hist2.bind_g[model]  # in erg
    bind_b = hist2.bind_b[model]
    bind_h = hist2.bind_h[model]
    print('bind_b = ', bind_b)

    m2 = msol * hist2.star_mass[model]
    a_i = recircul

    m2_core = msol * hist2.ce_core_mass[model]  # Chens ce_core_mass has H fraction 0.1 boundary condition
    # print(hist2.ce_core_mass[model], '---', hist2.he_core_mass[model])

    # envelope ejection efficiency:
    alpha = 0.5

    # calculate final separation in cm
    a_f = G * m2_core * M1rem / (2 * ((G * m2 * M1rem / (2 * a_i)) - (bind_b / alpha)))

    #q = m2_core/M1rem  # want Roche lobe radius of m2_core
    # TEST
    q = M1rem / m2_core  # see 2.37 in Daniels thesis

    # calculate Roche lobe radius in cm from Eggleton 1983
    R_L = (0.49 * (q ** (2 / 3)) * a_f) / (0.6 * (q ** (2 / 3)) + np.log(1 + q ** (1 / 3)))  # np.log = ln

    # what could be a good value that approximates the core radius? maybe event horizon of black hole companion?
    if a_f <= r_schw:
        print('The envelope cant be ejected before the two objects merge.')

    # compare final separation with Roche lobe radius of core

    elif a_f <= R_L:
        print('Final separation is smaller than Roche lobe radius of the naked core: merger')

    else:
        print('For alpha_CE = ', alpha)
        print('System survives common envelope.')
        print('Roche lobe of core in Rsun:', R_L / rsol)
        print('Final separation in Rsun:', a_f / rsol)
        P_orb = 2 * np.pi * np.sqrt((a_f ** 3) / (G * (M1rem + m2_core)))  # in seconds
        print('Final orbital period in days: ', P_orb / 86400)
        # TO DO: write final separation, Roche lobe and event horizon radius into an output file for the system

    # new calculation to get efficiency parameter from minimum final separation possible (R2_core = R_L)
    r2core = rsol * hist2.he_core_radius[model]

    # formula for Roche lobe radius used
    amin = r2core * (0.6 * (q ** (2 / 3)) + np.log(1 + q ** (1 / 3))) / (0.49 * (q ** (2 / 3)))

    Eorb = - (G * m2_core * M1rem)/(2*amin) + (G * m2 * M1rem)/(2 * a_i)
    # alphaCE = bind_b / (-(G * m2_core * M1rem / (2 * amin)) + (G * m2 * M1rem / (2 * a_i)))
    alphaCE = bind_b / Eorb

    print('efficiency parameter for minimal final binary separation (R2_core = R_RL): ', alphaCE)
    print('corresponding R2_core in Rsun: ', r2core / rsol)
    print('corresponding R2_core in cm: ', r2core)
    print('corresponding final binary separation in Rsun: ', amin / rsol)
    print('remaining core mass of secondary in Msun: ', m2_core / msol)
    print('remaining core mass of secondary in g: ', m2_core)
    print('envelope mass ejected during CE in Msun: ', (m2 - m2_core) / msol)


    # TO DO: any other parameters I need to start a new MESA run
    # C core already building in He star
    c_core_mass = hist2.c_core_mass[model]
    print('C core mass building in He star at point of CE: ', c_core_mass)

    # get orbital period from final binary separation
    P_orb = 2 * np.pi * np.sqrt((amin ** 3) / (G * (M1rem + m2_core)))  # in seconds
    print('Final orbital period in days for minimal separation: ', P_orb / 86400)

    print('system age where CE is happening in Myr: ', hist2.star_age[model] / 1e6)

    # TO DO: print neat file with the important values for the next MESA run