#
# Script for calculating Psi and Delta (Ellipsometry)
# or the Power of Transmission
# of variable amounts of periods
# with variable thickness
#
# More Materials may be added later.
#######################

import tmm
from numpy import pi, inf, arange
import numpy as np
import pandas as pd
import get_r_index
import matplotlib.pyplot as plt
# import calcback

degree = pi / 180
LAMBDAS = arange(400, 800+1, 1)  # Wavelength, minimum, maximum, number of steps
ANGLE = 70*degree

class layer:
    def __init__(self, material, d):
        self.material = material
        self.d = d
        self.calc_refractive_index(material)

    @classmethod
    def calc_refractive_index(cls, material):
        ''' gets the complex refractive index of the material for the desired wavelenghts from the get_r_index package'''
        # mat = get_r_index.material(material, LAMBDAS)
        # cls.cplx = np.around(mat.refractive_index, 3) #rounded to 3 decimals for faster computing

class period:
    def __init__(self, layer1, layer2):
        self.lay1 = layer1
        self.lay2 = layer2
        self.d = self.lay1.d + self.lay2.d


class stack:
    def __init__(self, period, num):
        self.period = period
        self.num = num # number of periods

    def get_d1(self):
        return self.period.lay1.d

    def get_d2(self):
        return self.period.lay2.d

# class sim:
#     def __init__(self, stack, lambdas=LAMBDAS):
#         self.stack = stack
#         # self.mode = mode # R, T, ellips
#         self.lambdas = lambdas
#         self.make_d_list(self.stack)

    # def e(self):
    #     self.ellips(self.stack)

def make_d_list(stack):
    d_list = [inf]
    for layer in range(0,stack.num):
        d_list.append(stack.get_d1())
        d_list.append(stack.get_d2())
    d_list.append(inf)
    return d_list

def T(stack):
    T_list = []
    for lambda_vac in LAMBDAS:
        n_list = [1]
        for layer in range(0,stack.num):
            n1 =  round(get_r_index.get_cplx(stack.period.lay1.material, lambda_vac),4)
            n_list.append(n1)
            n2 = round(get_r_index.get_cplx(stack.period.lay2.material, lambda_vac),4)
            n_list.append(n2)
        n_list.append(round(get_r_index.get_cplx('c-Si', lambda_vac), 4))
        T_list.append(tmm.coh_tmm('s', n_list, make_d_list(stack), ANGLE, lambda_vac)['T'])
    return T_list 

def R(stack):
    R_list = []
    for lambda_vac in LAMBDAS:
        n_list = [1]
        for layer in range(0,stack.num):
            n1 =  round(get_r_index.get_cplx(stack.period.lay1.material, lambda_vac),4)
            n_list.append(n1)
            n2 = round(get_r_index.get_cplx(stack.period.lay2.material, lambda_vac),4)
            n_list.append(n2)
        n_list.append(round(get_r_index.get_cplx('c-Si', lambda_vac), 4))
        R_list.append(tmm.coh_tmm('s', n_list, make_d_list(stack), ANGLE, lambda_vac)['R'])
    return R_list 

def ellips (stack):
    Psis = []
    Deltas = []
    for lambda_vac in LAMBDAS:
        n_list = [1]
        for layer in range(0,stack.num):
            n1 =  round(get_r_index.get_cplx(stack.period.lay1.material, lambda_vac),4)
            n_list.append(n1)
            n2 = round(get_r_index.get_cplx(stack.period.lay2.material, lambda_vac),4)
            n_list.append(n2)
        n_list.append(round(get_r_index.get_cplx('c-Si', lambda_vac), 4))
        e_data = tmm.ellips(n_list, make_d_list(stack), ANGLE, lambda_vac)
        Psis.append(e_data['psi']/degree)
        Deltas.append(e_data['Delta']/degree)
    output = np.asarray([Psis, Deltas])
    return output

def plot(dataframe, mode):
    df = dataframe
    if mode ==  1 or mode == 3:
        f1 = plt.figure(1)
        plt.ylabel('Psi[°]')
        # ax1.plot(df['lambda'], df['Psi'], label = 'Psi', color = 'tab:blue')
        f1 = df['Psi'].plot(color='b')
        f2 = plt.figure(2)
        f2 = df['Delta'].plot(color='r')
        plt.ylabel('Delta[°]')
        # ax2.plot(df['lamba'], df['Delta'], label = 'Delta', color = 'tab:red')
        plt.xlabel('Wavelength [nm]')
    if mode == 2 or mode == 3:
        plt.figure()
        fig2 = df['T'].plot()
        fig2.set_xlabel('Wavelength [nm]')
        fig2.set_ylabel('Power of Transmission')


def main():
    # Al = layer('Al', 0)
    # ITO = layer('ITO-RTA-II', 100)
    # aperiod = period(ITO, Al)
    ####
    # Ag = layer('Ag', 0)
    SiO2 = layer('SiO_2-II', 100)
    nothing = layer('SiO_2-II', 10) #we need 2 layers for a period, but thickness is 0
    aperiod = period(SiO2, nothing)
    thestack = stack(aperiod, 1)
    deltapsi = ellips(thestack)
    n_S = []
    for lambdavac in LAMBDAS:
        n_S.append(round(get_r_index.get_cplx('c-Si', lambdavac),4))
    df = pd.DataFrame({'lambda':LAMBDAS, 'Psi': deltapsi[0], 'Delta': deltapsi[1], 'n_S':np.real(n_S), 'k_S':np.imag(n_S)})
    df.set_index('lambda', inplace=True)
    plot(df, 1)
    df.to_csv('110nmSiO2.csv', encoding='utf-8')
    plt.show()


if __name__ == '__main__':
    main()
