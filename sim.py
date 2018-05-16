import tmm
# import tmm.examples
from numpy import pi, linspace, inf  # , array
# from scipy.interpolate import interpld
import matplotlib.pyplot as plt

degree = pi / 180
period = [10, 10]  # ammount, thickness (nm)
layer1 = 5  # standard

# Materials:
# Name = [n+jk, thickness in a period]
ITO = [1.635 + 0.0102j, layer1]
Al = [1.523 + 15.114j, (period[1] - ITO[1])]

lambda_list = linspace(300, 3000,
                       10000)  # Wavelength, minimum, maximum, number of steps
d_list = [inf]  # starts with air at infinity
for periods in range(0, period[0]):
    d_list.append(ITO[1])
    d_list.append(Al[1])
d_list.append(inf)  # end with wafer(?) at infinity

n_list = [1]  # starts with air
for periods in range(0, period[0]):
    n_list.append(ITO[0])
    n_list.append(Al[0])
n_list.append(3.49 + 0)  # Si (wafer)


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
    title = "{} x {} nm periods with {} nm ITO".format(period[0], period[1],
                                                       layer1)
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


def parameters():
    while True:
        try:
            period[0] = int(input('Number of periods: '))
            period[1] = int(input('Thickness of one period (nm): '))
            layer1 = int(input('Thickness of the ITO layer (nm): '))
            mode = input('Transmission (t) or Ellipsometry (e)? ')
            break
        except:
            print('Thats not a valid option')
    if mode == 't':
        transmission()
    elif mode == 'e':
        PsiDelta()
    else:
        print('Not supported, programm closed!')


if __name__ == "__main__":
    parameters()
    plt.show()
    print("Done")
