# supernova.py

import os
import numpy as np
import mesa_reader as mr
import matplotlib.pyplot as plt

path = os.getcwd()


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
# M1 = 64.0284 * msol  # ComBinE
# M2 = 38.1828 * msol  # ComBinE
M0 = M1 + M2  # total pre-supernova mass

# get remnant mass
McoreCO = msol * hist1.c_core_mass[model]  # c_core_mass
MlayerHe = msol * hist1.he_rich_layer_mass[model]  # he_rich_layer


M1rem = 0.8 * (McoreCO + 0.8 * MlayerHe)  # ComBinE recipe from line 1239 in evolution.cpp, see Kruckow 2018 2.2.6
print(M1rem)
# M1rem = 23.3292 * msol  # ComBinE
M1remsol = M1rem / msol

DeltaM = M1 - M1rem  # mass ejected from SN
a0 = rsol * histbi.binary_separation[model]  # pre-SN binary separation (binary_separation)

e0 = 0.0
r = a0 * (1 + (e0 ** 2) / 2)
w = 20.0 * (10 ** 5)  # kick velocity (from ComBinE screen output) in cm/s
theta = 90  # kick angle (from ComBinE screen output) in degrees
# np.cos needs input in rad
theta_rad = theta * np.pi / 180.0
phi = 90  # from ComBinE output
phi_rad = phi*np.pi/180.0




# test whether system becomes unbound by SN, according to J.G. Hills (1983)
# if more than half of system mass is lost by SN
if DeltaM / M0 >= (1 + (e0 ** 2)) / 2:
    print('System is destroyed by the SN')

else:
    vrel = (G * M0 / r) ** (1 / 2)


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



    Eorb = - (G * M1rem * M2 / r) + (mu * (vrel ** 2) / 2)


    e_final = np.sqrt(1 + (2 * Eorb * (Lorb ** 2)) / (mu * (G ** 2) * (M1rem ** 2) * (M2 ** 2)))

    # re-circularize orbit immediately (not realistic, but result should not be significantly different)
     recircul = (1 - e_final**2)*a
     recircul_rsol = (1 - e_final**2)*a_rsol


    # calculate Schwarzschild radius
    r_schw = 2 * G * M1rem / (c**2)


    # radiustest.py


    index = 0  # indexing of mesa reader lists
    lastmodel = hist2.model_number[-1]  # last model number, where for loop ends
    end = lastmodel - 1  # last index (starts with 0)
    print('last model of secondary:', lastmodel)

    # import radius and mass lists
    R2 = hist2.radius  # Rsun
    M2 = hist2.star_mass
    age = hist2.star_age
    # a0 is now input
    #######################################
    print('minimum post-SN orbital period in days:')
    P0_min = float(input())  # compare with Rsun value
    print('maximum post-SN orbital period in days:')
    P0_max = float(input())  # compare with Rsun value
    print('step size (int): ')
    step = int(input())
    num_steps = int((P0_max-P0_min)/step)
    print('#################################################################')

    Porb = []
    alpha_min = []

    for t in range(0, num_steps+1):
        P0 = P0_min + t*step
        Porb.append(P0)
        print('orbital period in days: ', P0)
        P0_yr = P0/365
        a0 = (G_sol * (M1remsol+M2[end]) * (P0_yr**2) / (4 * np.pi ** 2)) ** (1 / 3)
        print('post-SN binary separation in Rsun: ', a0)
        #######################################
        # make list with Roche lobe radii
        R_L = np.zeros(lastmodel)

        for i in range(index, end): # start at 0 for Roche lobe calculation
             # Roche lobe equation
            q = M2[i] / M1remsol
            # calculate Roche lobe radius in Rsun from Eggleton 1983
            RL = 0.49 * q ** (2 / 3) * a0 / (0.6 * q ** (2 / 3) + np.log(1 + q ** (1 / 3)))
            R_L[i] = RL
            if R2[i] >= RL:
                RL_last = RL
                crossover = i
                break

        # filling Roche lobe radius array with last calculated value for better indication
        # mention below figure
        for j in range(crossover, lastmodel):
            R_L[j] = RL_last

    # testing for convective mixing and mass ratio
    # largest region by mass is core

    # mass ratio
        ratio = M2[crossover] / M1remsol
        print('mass ratio at RLO: ', ratio)
        print('model number: ', crossover)

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
            print('System survives common envelope.')

            P_orb = 2 * np.pi * np.sqrt((a_f ** 3) / (G * (M1rem + m2_core)))  # in seconds

     # new calculation to get efficiency parameter from minimum final separation possible (R2_core = R_L)
        r2core = rsol * hist2.he_core_radius[model]

    # formula for Roche lobe radius used
        amin = r2core * (0.6 * (q ** (2 / 3)) + np.log(1 + q ** (1 / 3))) / (0.49 * (q ** (2 / 3)))

        Eorb = - (G * m2_core * M1rem) / (2 * amin) + (G * m2 * M1rem) / (2 * a_i)
        # alphaCE = bind_b / (-(G * m2_core * M1rem / (2 * amin)) + (G * m2 * M1rem / (2 * a_i)))
        alphaCE = bind_b / Eorb

        print('efficiency parameter for minimal final binary separation (R2_core = R_RL): ', alphaCE)
        alpha_min.append(alphaCE)

        print('#################################################################')



plot = plt.figure(dpi=600)

plt.xlabel('$P_{\mathrm{orb, pre-CE}}\ /\ \mathrm{days}$')
plt.ylabel('$\\alpha_{\mathrm{CE, min}}$')
plt.tick_params(axis='both', direction='in', which='both', top=True, right=True)
plt.ylim(0, 0.6)
plt.axes().set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])

plt.plot(Porb, alpha_min, marker='', ls='dashed', c='green')




#plt.title('Radii of the stellar components over time with a comparison between MESA and ComBinE results', wrap=True)

plt.savefig('alpha_min')
plt.savefig('alpha_min.pdf')

plt.show()



