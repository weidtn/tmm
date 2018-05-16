##############################################################################
#
# Script for calculating Psi and Delta (Ellipsometry)
# or the Power of Transmission
# of a single layer of ITO on top of SiO2
#
###################################

import tmm
import tmm.examples
from numpy import pi, linspace, inf  # , array
import matplotlib.pyplot as plt

degree = pi / 180

ITO = [1.635 + 0.0102j, 100]
SiO2 = [1.44, 380]
Si = [3.67]

lambda_list = linspace(300, 2000,
                       10000)  # Wavelength, minimum, maximum, number of steps
d_list = [inf, 100, 300, inf]
n_list = [1, ITO[0], SiO2[0], Si[0]]

def PsiDelta():
    psis = []
    Deltas = []
    for lambda_vac in lambda_list:
        e_data = tmm.tmm_core.ellips(n_list, d_list, 70 * degree, lambda_vac)
        psis.append(e_data['psi'] / degree)  # angle in degrees
        Deltas.append(e_data['Delta'] / degree)  # angle in degrees
    plt.figure()
    plt.plot(lambda_list, psis, lambda_list, Deltas)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Psi and Delta')
    # plt.title('test1')


def transmission():

    T_list = []
    for lambda_vac in lambda_list:
        T_list.append(
            tmm.tmm_core.coh_tmm('s', n_list, d_list, 0, lambda_vac)['T'])
    plt.figure()
    plt.plot(lambda_list, T_list)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Fraction of power Transmitted')

# transmission()
PsiDelta()
plt.show()
print("Done")
