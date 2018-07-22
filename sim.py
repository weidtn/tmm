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
from numpy import pi, inf, linspace
import get_r_index
import matplotlib.pyplot as plt

degree = pi / 180
lambda_list = linspace(450, 3000,
                       10000)  # Wavelength, minimum, maximum, number of steps


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

    aperiod = period(d, layer1)

    # Setting up the list with thicknesses
    d_list = [inf]  # starts with air at infinity
    for periods in range(0, num):
        d_list.append(aperiod.get_layer1())
        if aperiod.get_layer2() != 0:
            d_list.append(aperiod.get_layer2())
    d_list.append(inf)  # end with wafer(?) at infinity

    # # Starting the calculations:
    if mode == 't':
        T_list = []
        for lambda_vac in lambda_list:
            # Setting up the list with indices of refraction
            n_list = [1]  # starts with air
            for periods in range(0, num):
                n_list.append(get_r_index.get_cplx('ITO', lambda_vac))
                if aperiod.get_layer2() not 0:
                    n_list.append(get_r_index.get_cplx('Al', lambda_vac))
                n_list.append(3.49 + 0)  # Si (wafer)
            T_list.append(
                tmm.coh_tmm('s', n_list, d_list, 70 * degree, lambda_vac)['T'])
        plt.figure()
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Fraction of power Transmitted')
        plt.plot(lambda_list, T_list)

    elif mode == 'e':
        psis = []
        Deltas = []
        for lambda_vac in lambda_list:
            # Setting up the list with indices of refraction
            n_list = [1]  # starts with air
            for periods in range(0, num):
                n_list.append(get_r_index.get_cplx('ITO', lambda_vac))
                if aperiod.get_layer2() not 0:
                    n_list.append(get_r_index.get_cplx('Al', lambda_vac))
            n_list.append(3.49 + 0)  # Si (wafer)
            e_data = tmm.tmm_core.ellips(n_list, d_list, 70 * degree,
                                         lambda_vac)
            psis.append(e_data['psi'] / degree)  # angle in degrees
            Deltas.append(e_data['Delta'] / degree)  # angle in degrees
        plt.figure()
        plt.plot(lambda_list, psis, lambda_list, Deltas)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Psi and Delta')
        title = "{} x {} nm periods with {} nm ITO".format(
            int((len(d_list) - 2) / 2), d_list[1] + d_list[2], d_list[1])
        plt.title(title)

    elif mode == 'debug':
        print('debug')
        print(n_list)
    else:
        print('Not supported, programm closed!')


if __name__ == "__main__":
    main()
    plt.show()
    print("Done")
