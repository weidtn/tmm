######################
#
# Script for calculating Psi and Delta (Ellipsometry)
# or the Power of Transmission
# of variable ammount of periods
# made out of alternating layers of ITO and Al
# with variable thickness
#
# More Materials may be added later.
#######################

import tmm
from numpy import pi, linspace, inf
import matplotlib.pyplot as plt

degree = pi / 180


# Materials:
# Name = [n+jk, thickness in a period]
class material:
    def __init__(self, n):
        # n is the complex index of refraction
        # given as 0.000+0.000j
        self.n = n

    def get_n(self):
        return self.n


class period:
    def __init__(self, d, layer1):
        self.d = d  # Thickness of one Period
        self.layer1 = layer1  # Thickness of the first layer
        self.layer2 = (d - layer1)

    def get_layer1(self):
        return self.layer1

    def get_layer2(self):
        return self.layer2

    def get_d(self):
        return self.d


def PsiDelta(n_list, d_list, lambda_list):
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
    title = "{} x {} nm periods with {} nm ITO".format(
        int((len(d_list) - 2) / 2), d_list[1] + d_list[2], d_list[1])
    plt.title(title)


def transmission():

    T_list = []
    for lambda_vac in lambda_list:
        T_list.append(
            tmm.tmm_core.coh_tmm('s', n_list, d_list, 0, lambda_vac)['T'])
    plt.figure()
    plt.plot(lambda_list, T_list)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Fraction of power Transmitted')


# ---------- Main function: ------------------


def main():
    while True:
        try:
            num = int(input('Number of periods: '))
            d = int(input('Thickness of one period (nm): '))
            layer1 = int(
                input('Thickness of the first layer in a period (nm): '))
            mode = input('Transmission (t) or Ellipsometry (e)? ')
            break
        except:
            print('Thats not a valid option!')

    lambda_list = linspace(
        450, 3000, 10000)  # Wavelength, minimum, maximum, number of steps

    aperiod = period(d, layer1)
    ITO = material(1.64 + 0.01j)
    Al = material(1.52 + 15.1j)
    # Setting up the list with thicknesses
    d_list = [inf]  # starts with air at infinity
    for periods in range(0, num):
        d_list.append(aperiod.get_layer1())
        d_list.append(aperiod.get_layer2())
    d_list.append(inf)  # end with wafer(?) at infinity

    # Setting up the list with indices of refraction
    n_list = [1]  # starts with air
    for periods in range(0, num):
        n_list.append(ITO.get_n())
        n_list.append(Al.get_n())
    n_list.append(3.49 + 0)  # Si (wafer)

    # Starting the calculations:
    if mode == 't':
        transmission()
    elif mode == 'e':
        PsiDelta(n_list, d_list, lambda_list)
    elif mode == 'debug':
        print('debug')
        print(type(aperiod.get_layer1()))
    else:
        print('Not supported, programm closed!')


if __name__ == "__main__":
    main()
    plt.show()
    print("Done")
