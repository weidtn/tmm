import tmm
import tmm.examples
from numpy import pi, linspace, inf  # , array
# from scipy.interpolate import interpld
import matplotlib.pyplot as plt

# Single-layer film with complicated wavelength-dependent n :
# tmm.examples.sample2()

degree = pi / 180
period = [100, 10]  # ammount, thickness (nm)
layer1 = input('Thickness of the first ITO layer: ')

# Materials:
# Name = [n+jk, thickness in a period]
ITO = [1.635 + 0.0102j, int(layer1)]
Al = [1.523 + 15.114j, (period[1] - ITO[1])]

# for recreating the matlab-script for InP-Air-DBR
# inp = [122.355]
# air = 387.5
# d_list = [
#     inf, inp, air, inp, air, inp, 815, inp, air, inp, air, inp, air, inf
# ]  # Thicknesses in nm
# n_list = [1, 3.167, 1, 3.167, 1, 3.167, 1, 3.167, 1, 3.167, 1, 3.167, 1, 3.167]

lambda_list = linspace(1000, 3000,
                       10000)  # Wavelength, minimum, maximum, number of steps
d_list = [inf]  # starts with air at infinity
for periods in range(0, period[0]):
    d_list.append(ITO[1])
    d_list.append(Al[1])
d_list.append(inf) # end with wafer? at infinity

n_list = [1] # starts with air
for periods in range(0,period[0]):
    n_list.append(ITO[0])
    n_list.append(Al[0])
n_list.append(3.49+0)  # Si (wafer)


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
