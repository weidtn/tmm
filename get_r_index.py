"""
Script to calculate refractive indices for a range of wavelengths

-----------------------------------------
INPUT variables:

mat_name = Material name (string)
wl = Wavelength as numpy array
imSign = n+jk: '+'
         n-jk: '-'
         (optional, standard is +)

-------------------------------
Your material can OUTPUT these variables (numpy arrays):

refractive_index
n
k
epsilon
epsilon_real
epsilon_imag
wl

---------------------------------
Available materials for these models:

Drude-Lorentz:  Au, ITO, ITO-RTA, Ag, Al, Cu, Cr, Ni, W, Ti, Pt
Tauc-Lorentz: SiO_2, SiO_2-II, Si_3N_4, ZrO_2, ZrO_2-II, a-Si
Sellmeier: no materials yet, but function should work
Leng-Lorentz: c-Si // TODO something is not right here
-----------------
Example usage in your scripts:
import get_r_index
Au = get_r_index.material(material='Au', wl=numpy.arange(100,1000,1), imSign='+')
print(Au.refractive_index)

OR:
Au = get_r_index.get_cplx(material='Au', wl=1000)      # Gets the complex refractive index for ONE wavelength
print(Au)
-------------------
You can also run this file to quickly plot n and k for


"""
import numpy as np
import matplotlib.pyplot as plt

SUPPORTED_MATERIALS = [
    'Au', 'ITO', 'ITO-RTA', 'Ag', 'Al', 'Cu', 'Cr', 'Ni', 'W', 'Ti', 'Pt',
    'SiO_2', 'SiO_2-II', 'Si_3N_4', 'ZrO_2', 'ZrO_2-II', 'a-Si'
]


class material():
    def __init__(self, mat_name, wl, imSign='+'):

        self.mat_name = mat_name
        self.wl = wl
        self.imSign = imSign

        if self.mat_name == 'SiO_2':
            # fitting coefficients (IBD)
            # B1 = -27.424
            # B2 = 1.3995
            # B3 = 27.23
            # C1 = 5643.3
            # C2 = 10805
            # C3 = 5407.7
            # vali = [400, 2000]
            # model = 'Sellmeier'
            #
            #        # applied material model:  Sellmeier  (only VIS + NIR)
            #        #[cplx, reP, imP] = n_Sellmeier(wl, B1, B2, B3, C1, C2, C3, RIorEPS, imSign, vali)

            #        # fitting coefficients (IBD Taimoor)
            #        C0 = 100
            #        C1 = 10000000
            #        N0 = 1.482
            #        N1 = 54.9
            #        N2 = -29
            #        K0 = 0
            #        K1 = 0
            #        K2 = 0
            #        vali = [400 1000]
            #
            #        # applied material model: Cauchy
            #        #[cplx, reP, imP] = n_Cauchy(wl, C0, C1, N0, N1, N2, K0, K1, K2, RIorEPS, imSign, vali)

            # # fitting coefficients  (IBD Hans)
            EgeV = 7.5554
            AeV = 106.67
            E0eV = 9.9221
            CeV = 0.3209
            eps_inf = 1.57696
            vali = [100, 10000]
            self.model = 'Tauc_Lorentz'
        # applied material model: Tauc-Lorentz
        #[cplx, reP, imP] = n_Tauc_Lorentz(wl, EgeV, E0eV, AeV, CeV, eps_inf, RIorEPS, imSign, vali)

        elif self.mat_name == 'SiO_2-II':

            # fitting coefficients  (IBD Taimoor)
            EgeV = 7.56
            AeV = 77.4
            E0eV = 9.9639
            CeV = 0.7802
            eps_inf = 1.72683
            vali = [100, 10000]
            self.model = 'Tauc_Lorentz'

        # applied material model: Tauc-Lorentz
        #[cplx, reP, imP] = n_Tauc_Lorentz(wl, EgeV, E0eV, AeV, CeV, eps_inf,                                          RIorEPS, imSign, vali)

        elif self.mat_name == 'Si_3N_4':

            # fitting coefficients (Palik, Handbook of Optical Constants of Solids, AP)
            EgeV = 4.5
            AeV = 59.2
            E0eV = 6.78
            CeV = 0.49
            eps_inf = 3.1
            vali = [100, 10000]
            self.model = 'Tauc_Lorentz'

        # applied material model: Tauc-Lorentz
        #[cplx, reP, imP] = n_Tauc_Lorentz(wl, EgeV, E0eV, AeV, CeV, eps_inf,                                          RIorEPS, imSign, vali)

        elif self.mat_name == 'Nb_2O_5':

            #         # fitting coefficients (IBD)
            #         B1 = 1.5546
            #         B2 = 1.559
            #         B3 = 0.86685
            #         C1 = 54525
            #         C2 = 54527
            #         C3 = -2.3879e5
            #         vali = [400 2000]
            #
            #         # applied material model:  Sellmeier
            #          #[cplx, reP, imP] = n_Sellmeier(wl, B1, B2, B3, C1, C2, C3, RIorEPS, imSign, vali)

            # fitting coefficients (IBD Taimoor)
            C0 = 100
            C1 = 10000000
            N0 = 2.225
            N1 = 180.6
            N2 = 377
            K0 = 0
            K1 = 0
            K2 = 0
            vali = [400, 1000]

        # applied material model: Cauchy
        #[cplx, reP, imP] = n_Cauchy(wl, C0, C1, N0, N1, N2, K0, K1, K2,                                    RIorEPS, imSign, vali)
        elif self.mat_name == 'NCD':

            # fitting coefficients (CVD Cyril)
            C0 = 100
            C1 = 10000000
            N0 = 2.397
            N1 = 72.2
            N2 = 56.6
            K0 = 0
            K1 = 34.303
            K2 = 0.002
            vali = [400, 1000]

        # applied material model: Cauchy
        #[cplx, reP, imP] = n_Cauchy(wl, C0, C1, N0, N1, N2, K0, K1, K2,                                    RIorEPS, imSign, vali)

        elif self.mat_name == 'Al_2O_3':

            #           # fitting coefficients (IBD)
            #           B1 = 5.17502
            #           B2 = -3.40698
            #           B3 = -79.6032
            #           C1 = 1.10646e4
            #           C2 = 1.10363e4
            #           C3 = -7.51505e9
            #           vali = [400 2000]
            #
            #           # applied material model:  Sellmeier
            #           #[cplx, reP, imP] = n_Sellmeier(wl, B1, B2, B3, C1, C2, C3, RIorEPS, imSign, vali)

            # fitting coefficients (Horiba, New Amorphous Dispersion Formula, Technical Note (2006))
            n_inf = 1.56
            w_g = 9.85
            amp = [2.46]
            frq = [10.4]
            dmp = [0.44]
            vali = [250, 1700]

        # applied material model: Horiba New Amorphous
        #[cplx, reP, imP] = n_Horiba_New_Amorphous(wl, n_inf, w_g, amp, frq,                                                  dmp, RIorEPS, imSign, vali)

        elif self.mat_name == 'TiO_2':

            #           # fitting coefficients (IBD)
            #           B1 = -499.79
            #           B2 = 2.113
            #           B3 = 501.73
            #           C1 = -1.3658e5
            #           C2 = 82237
            #           C3 = -1.357e5
            #           vali = [400 2000]
            #
            #           # applied material model: Sellmeier
            #           [n_cplx, n_mat, k_mat] = n_Sellmeier(wl, B1, B2, B3, C1, C2, C3, RIorEPS, imSign, vali)

            # fitting coefficients (Horiba, New Amorphous Dispersion Formula, Technical Note (2006))
            n_inf = 2.096
            w_g = 2.952
            amp = [0.278]
            frq = [4.282]
            dmp = [0.652]
            vali = [300, 1700]

        # applied material model: Horiba New Amorphous
        #[cplx, reP, imP] = n_Horiba_New_Amorphous(wl, n_inf, w_g, amp, frq,dmp, RIorEPS, imSign, vali)

        elif self.mat_name == 'ZrO_2':

            #         # fitting coefficients (IBD)
            #         B1 = -7.3043
            #         B2 = -1.3686
            #         B3 = 12.177
            #         C1 = -2296
            #         C2 = 20571
            #         C3 = 10609
            #         vali = [400 2000]
            #
            #         # applied material model: Sellmeier
            #         [n_cplx, n_mat, k_mat] = n_Sellmeier(wl, B1, B2, B3, C1, C2, C3, RIorEPS, imSign, vali)

            # fitting coefficients (IBD Hans)
            EgeV = 3.6053
            AeV = 111.17
            E0eV = 4.2382
            CeV = 2.7505
            eps_inf = 3.32592
            vali = [100, 10000]
            self.model = 'Tauc_Lorentz'

        # applied material model: Tauc-Lorentz
        #[cplx, reP, imP] = n_Tauc_Lorentz(wl, EgeV, E0eV, AeV, CeV, eps_inf,                                          RIorEPS, imSign, vali)

        elif self.mat_name == 'ZrO_2-II':

            # fitting coefficients (IBD Taimoor)
            EgeV = 4.7926
            AeV = 359.08
            E0eV = 5.1421
            CeV = 2.7505
            eps_inf = 2.78023
            vali = [100, 10000]
            self.model = 'Tauc_Lorentz'

        # applied material model: Tauc-Lorentz
        #[cplx, reP, imP] = n_Tauc_Lorentz(wl, EgeV, E0eV, AeV, CeV, eps_inf,                                          RIorEPS, imSign, vali)

    # elif self.mat_name == 'SodaLime':

    #     # fitting coefficients (based on: Synowicki et al., Thin Solid Films, 519/9, 2907-2913 (2011))
    #     order = 9          # polynome order
    #     cutoff = [0 500]       # transition between exponential and polynome function for n and k
    #     en_par = [0 1 0]        # exponetial function parameters (amplitude, decay time, offset) for n
    #     ek_par = [2340.939752268 17.420239546489 1e-07]        # exponetial function parameters (amplitude, decay time, offset) for k
    #     pn_amp = [1.8911444029173 -0.002585073761333 7.9982090232581e-06 ...
    #               -1.4322974036864e-08 1.6015031221774e-11 -1.1512147623537e-14 ...
    #               5.3139948976884e-18 -1.5217369316064e-21 2.4605040747393e-25 ...
    #               -1.716994764055e-29]        # polynome amplitudes, length: order+1, for n
    #     pk_amp = [-0.00032417187842718 2.594721462676e-06 -8.7227626441111e-09 ...
    #               1.6125995711296e-11 -1.8080697492943e-14 1.2825310420403e-17 ...
    #               -5.7954637604999e-21 1.6191916440455e-24 -2.5523663790964e-28 ...
    #               1.7377319738944e-32]        # polynome amplitudes, length: order, for k
    #     vali = [300 2500]

    #     # applied material model: polynome + exponential fit
    #     #[cplx, reP, imP] = n_polynome_exp(wl, order, cutoff, en_par, ek_par, pn_amp, pk_amp, RIorEPS, imSign, vali)

        elif self.mat_name == 'ITO':

            # fitting coefficients (IBD, Wetterau)
            w_p = 1.486  # plasma freq
            amp = [1, 43.97]  # oscillator strength
            gam = [0.148, 0.945]  # damping
            frs = [0, 4.655]  # resonance freq
            vali = [100, 10000]
            self.model = 'Drude_Lorentz'

        elif self.mat_name == 'ITO-RTA':

            # fitting coefficients (IBD, Wetterau)
            w_p = 1.944  # plasma freq
            amp = [1, 32.248]  # oscillator strength
            gam = [0.206, 0.274]  # damping
            frs = [0, 5.136]  # resonance freq
            vali = [100, 10000]
            self.model = 'Drude_Lorentz'

        elif self.mat_name == 'ITO-II':

            # fitting coefficients (IBD, Künne)
            w_p = 1.415  # plasma freq
            amp = [1, 47.4141]  # oscillator strength
            gam = [0.062, 0.487]  # damping
            frs = [0, 6.156]  # resonance freq
            vali = [100, 10000]
            self.model = 'Drude_Lorentz'

        elif self.mat_name == 'ITO-RTA-II':

            # fitting coefficients (IBD, Künne)
            w_p = 1.932  # plasma freq
            amp = [1, 4.402]  # oscillator strength
            gam = [0.025, 0.053]  # damping
            frs = [0, 5.011]  # resonance freq
            vali = [100, 10000]
            self.model = 'Drude_Lorentz'

        elif self.mat_name == 'Au':

            # fitting coefficients (Raki? et al., Appl. Opt. 37, 5271-5283 (1998))
            w_p = 9.03  # plasma freq
            amp = [0.76, 0.024, 0.01, 0.071, 0.601,
                   4.384]  # oscillator strength
            gam = [0.053, 0.241, 0.345, 0.87, 2.494, 2.214]  # damping
            frs = [0, 0.415, 0.83, 2.969, 4.304, 13.32]  # resonance freq
            vali = [100, 10000]
            self.model = 'Drude_Lorentz'

        elif self.mat_name == 'Ag':

            # fitting coefficients (Raki? et al., Appl. Opt. 37, 5271-5283 (1998))
            w_p = 9.01  # plasma freq
            amp = [0.845, 0.065, 0.124, 0.011, 0.84,
                   5.646]  # oscillator strength
            gam = [0.048, 3.886, 0.452, 0.065, 0.916, 2.419]  # damping
            frs = [0, 0.816, 4.481, 8.185, 9.083, 20.29]  # resonance freq
            vali = [100, 10000]
            self.model = 'Drude_Lorentz'

        elif self.mat_name == 'Al':

            # fitting coefficients (Raki? et al., Appl. Opt. 37, 5271-5283 (1998))
            w_p = 14.98  # plasma freq
            amp = [0.523, 0.227, 0.050, 0.166, 0.030]  # oscillator strength
            gam = [0.047, 0.333, 0.312, 1.351, 3.382]  # damping
            frs = [0.000, 0.162, 1.544, 1.808, 3.473]  # resonance freq
            vali = [100, 10000]
            self.model = 'Drude_Lorentz'

        elif self.mat_name == 'Cu':

            # fitting coefficients (Raki? et al., Appl. Opt. 37, 5271-5283 (1998))
            w_p = 10.83  # plasma freq
            amp = [0.575, 0.061, 0.104, 0.723, 0.638]  # oscillator strength
            gam = [0.030, 0.378, 1.056, 3.213, 4.305]  # damping
            frs = [0.000, 0.291, 2.957, 5.300, 11.18]  # resonance freq
            vali = [100, 10000]
            self.model = 'Drude_Lorentz'

        elif self.mat_name == 'Cr':

            # fitting coefficients (Raki? et al., Appl. Opt. 37, 5271-5283 (1998))
            w_p = 10.75  # plasma freq
            amp = [0.168, 0.151, 0.150, 1.149, 0.825]  # oscillator strength
            gam = [0.047, 3.175, 1.305, 2.676, 1.335]  # damping
            frs = [0.000, 0.121, 0.543, 1.970, 8.775]  # resonance freq
            vali = [100, 10000]
            self.model = 'Drude_Lorentz'

        elif self.mat_name == 'Ni':

            # fitting coefficients (Raki? et al., Appl. Opt. 37, 5271-5283 (1998))
            w_p = 15.92  # plasma freq
            amp = [0.096, 0.100, 0.135, 0.106, 0.729]  # oscillator strength
            gam = [0.048, 4.511, 1.334, 2.178, 6.292]  # damping
            frs = [0.000, 0.174, 0.582, 1.597, 6.089]  # resonance freq
            vali = [100, 10000]
            self.model = 'Drude_Lorentz'

        elif self.mat_name == 'W':

            # fitting coefficients (Raki? et al., Appl. Opt. 37, 5271-5283 (1998))
            w_p = 13.22  # plasma freq
            amp = [0.206, 0.054, 0.166, 0.706, 2.590]  # oscillator strength
            gam = [0.064, 0.530, 1.281, 3.332, 5.836]  # damping
            frs = [0.000, 1.004, 1.917, 3.580, 7.498]  # resonance freq
            vali = [100, 10000]
            self.model = 'Drude_Lorentz'

        elif self.mat_name == 'Ti':

            # fitting coefficients (Raki? et al., Appl. Opt. 37, 5271-5283 (1998))
            w_p = 7.29  # plasma freq
            amp = [0.148, 0.899, 0.393, 0.187, 0.001]  # oscillator strength
            gam = [0.082, 2.276, 2.518, 1.663, 1.762]  # damping
            frs = [0.000, 0.777, 1.545, 2.509, 1.943]  # resonance freq
            vali = [100, 10000]
            self.model = 'Drude_Lorentz'

        elif self.mat_name == 'Pt':

            # fitting coefficients (Raki? et al., Appl. Opt. 37, 5271-5283 (1998))
            w_p = 9.59  # plasma freq
            amp = [0.333, 0.191, 0.659, 0.547, 3.576]  # oscillator strength
            gam = [0.080, 0.517, 1.838, 3.668, 8.517]  # damping
            frs = [0.000, 0.780, 1.314, 3.141, 9.249]  # resonance freq
            vali = [100, 10000]
            self.model = 'Drude_Lorentz'

        elif self.mat_name == 'c-Si':
            # fitting coefficients (Sentech SENpro Software)
            C0 = [56.32, 240.9, 125.46, 16.66]
            Beta = [-0.4589, -0.411, 0.3307, 0.2816]
            Eg = [3.38, 3.6266, 4.2906, 5.3825]
            Gam = [0.115, 0.3079, 0.203, 0.241]
            Mu = [-0.8241, -0.3965, -0.9504, -0.9761]
            eps_inf = 0.397200
            m0 = -0.033029
            x0 = 1.012892
            k0 = 0.0241
            vali = [100, 10000]
            self.model = 'Leng_Lorentz'

        elif self.mat_name == 'a-Si':

            # fitting coefficients (Palik, Handbook of Optical Constants of Solids, AP)
            EgeV = 1.2
            AeV = 122
            E0eV = 3.45
            CeV = 2.54
            eps_inf = 1.15
            vali = [100, 10000]
            self.model = 'Tauc_Lorentz'

        # applied material model: Tauc-Lorentz
        ##[cplx, reP, imP] = n_Tauc_Lorentz(wl, EgeV, E0eV, AeV, CeV, eps_inf,                                          RIorEPS, imSign, vali)

        elif self.mat_name == 'PMMA':

            # fitting coefficients (Horiba, New Amorphous Dispersion Formula, Technical Note (2006))
            n_inf = 1.456
            w_g = 3.667
            amp = [0.13]
            frq = [4.212]
            dmp = [7.144]
            vali = [200, 2100]

        # applied material model: Horiba New Amorphous
        ##[cplx, reP, imP] = n_Horiba_New_Amorphous(wl, n_inf, w_g, amp, frq,                                                  dmp, RIorEPS, imSign, vali)

        elif self.mat_name == 'InP':

            # fitting coefficients (MOCVD)
            eps_f0 = 4.8
            fres = [1.24, 2.62]
            gamm = [0.001, 0.29]
            sigm = [0.33, 4.4]
            vali = [1000, 2000]

        # applied material model: Lorentz (meep), (only NIR > 1000 nm)
        ##[cplx, reP, imP] = n_meep_Lorentz(wl, eps_f0, fres, gamm, sigm,                                         RIorEPS, imSign, vali)

        else:
            print('Error, material name not supported! Try \'list\' ')

        if self.model == 'Drude_Lorentz':
            self.n_Drude_Lorentz(wl, w_p, amp, gam, frs, vali)
        elif self.model == 'Tauc_Lorentz':
            self.n_Tauc_Lorentz(wl, EgeV, AeV, E0eV, CeV, eps_inf, vali)
        elif self.model == 'Sellmeier':
            self.n_Sellmeier(wl, B1, B2, B3, C1, C2, C3, vali)
        elif self.model == 'Leng_Lorentz':
            self.n_Leng_Lorentz(wl, C0, Beta, Eg, Gam, Mu, eps_inf, m0, x0, k0,
                                vali)
        else:
            print('Not a valid model!')

#########################################################

# ========== dispersion model functions
# function ##[cplx, reP, imP] = n_Cauchy(lambda, C0, C1, N0, N1, N2, K0, K1, K2, RIorEPS, imSign, vali)

#           if ((min(lambda) < vali(1)) || (max(lambda) > vali(2)))
#               warning('Wavelength out ')
#               cplx = reP + 1i*imP;
#           else
#               warning('Unknown option for parameter sign');
#           end

    @classmethod
    def n_Sellmeier(cls, wl, B1, B2, B3, C1, C2, C3, vali):

        if wl[0] < vali[0] or wl[1] > vali[1]:
            print('Wavelenght out of model validity range!')
        lambdasquare = np.squared(wl)
        n_squared = 1 + np.divide(
            (B1 * lambdasquare), (lambdasquare - C1)) + np.divide(
                (B2 * lambdasquare), (lambdasquare - C2)) + np.divide(
                    (B3 * lambdasquare), (lambdasquare - C3))
        n_mat = np.sqrt(n_squared)
        #           n_mat = sqrt(1 + B1*lambda.^2 ./ (lambda.^2 - C1) + B2*lambda.^2 ./ (lambda.^2 - C2) + B3*lambda.^2 ./ (lambda.^2 - C3));
        k_mat = np.divide(1, (np.multiply(
            n_mat,
            (B1 * wl + np.divide(B2, wl) + np.divide(B3, np.power(wl, 3))))))
        #           k_mat = 1./(n_mat .* (B1*lambda + B2./lambda + B3./lambda.^3));
        cls.n = n_mat
        cls.k = k_mat
        cls.refractive_index = cls.n + cls.k * 1j
        cls.epsilon_real = np.square(n_mat) - np.square(k_mat)
        cls.epsilon_imag = 2 * np.multiply(n_mat, k_mat)

    @classmethod
    def n_Drude_Lorentz(cls, wl, w_p, amp, gam, frs, vali):

        # --- for pure Drude set only amp(1) non-zero
        # --- for pure Lorentz set only amp(1) to zero
        if wl[0] < vali[0] or wl[1] > vali[1]:
            print('Wavelenght out of model validity range!')

        # convert wavelength to eV
        EeV = 1240 / wl

        epsC = 1 - np.divide((amp[0] * w_p**2),
                             (EeV**2 + gam[0] * EeV * 1j))  # Drude term:

        for i in range(0, len(amp), 1):  #  Lorentz terms
            epsC += (amp[i] * w_p**2) / (
                frs[i]**2 - EeV**2 - gam[i] * EeV * 1j)

        n = np.sqrt(
            0.5 * (np.sqrt((epsC.real**2) + (epsC.imag**2)) + epsC.real))
        k = np.sqrt(
            0.5 * (np.sqrt((epsC.real**2) + (epsC.imag**2)) - epsC.real))
        cls.epsilon = epsC
        cls.epsilon_real = epsC.real
        cls.epsilon_imag = epsC.imag
        cls.n = n
        cls.k = k
        cls.refractive_index = n + k * 1j

    @classmethod
    def n_Leng_Lorentz(cls, wl, C0, Beta, Eg, Gam, Mu, eps_inf, m0, x0, k0,
                       vali):

        # Source: Leng et al., JVST A, 16:3, 1654-1657 (1998)
        if wl[0] < vali[0] or wl[1] > vali[1]:
            print('Wavelenght out of model validity range!')
        # convert lambda to eV
        EeV = 1240 / wl
        epsC = eps_inf + m0 * EeV**x0 + 1j * k0 

        for i in range(0, len(C0), 1):
            epsC = epsC + (C0[i] / EeV**2) * (
                np.exp(1j * Beta[i]) *
                (Eg[i] - EeV - 1j * Gam[i])**Mu[i] + np.exp(-1j * Beta[i]) *
                (Eg[i] + EeV + 1j * Gam[i])**Mu[i]) - 2 * np.real(
                    np.exp(-1j * Beta[i]) * (Eg[i] + 1j * Gam[i])**Mu[i]
                ) - 2 * 1j * Mu[i] * EeV * np.imag(
                    np.exp(-1j * Beta[i]) * (Eg[i] + 1j * Gam[i])**(Mu[i] - 1))

        n = np.sqrt(
            0.5 * (np.sqrt((epsC.real**2) + (epsC.imag**2)) + epsC.real))
        k = np.sqrt(
            0.5 * (np.sqrt((epsC.real**2) + (epsC.imag**2)) - epsC.real))
        cls.epsilon = epsC
        cls.epsilon_real = epsC.real
        cls.epsilon_imag = epsC.imag
        cls.n = n
        cls.k = k
        cls.refractive_index = n + k * 1j

        # function ##[cplx, reP, imP] = n_meep_Lorentz(lambda, eps_f0, fres, gamm, sigm, RIorEPS, imSign, vali)


#           if ((min(lambda) < vali(1)) || (max(lambda) > vali(2)))
#               warning('Wavelength out of model validity range!');
#           end

#           # convert lambda to normalized frequency (meep a=1)
#           frq = 1000./lambda;

#           epsC = eps_f0;

#           for idx = 1:1:length(fres)
#               epsC = epsC + ((fres(idx).^2 .* sigm(idx)) ./ (fres(idx).^2 - frq.^2 - 1i*gamm(idx) .* frq));
#           end

    @classmethod
    def n_Tauc_Lorentz(cls, wl, EgeV, AeV, E0eV, CeV, eps_inf, vali):

        # Source: Jellison & Modine, APL 69:3, 371-373 (1996) and
        # Erratum: Jellison & Modine, APL 69:14, 2137 (1996)

        if wl[0] < vali[0] or wl[1] > vali[1]:
            print('Wavelenght out of model validity range!')

        Eg = EgeV
        E0 = E0eV
        A = AeV
        C = CeV
        epsR = eps_inf
        epsI = 0
        E = (1240 / wl)

        # Substitutions:
        t1 = (Eg**2 - E0**2) * E**2
        t2 = (Eg**2 * C**2)
        t3 = E0**2 * (E0**2 + 3 * Eg**2)
        a_ln = t1 + t2 - t3

        t1 = (E**2 - E0**2) * (E0**2 + Eg**2)
        t2 = Eg**2 * C**2
        a_atan = t1 + t2

        alpha = np.sqrt(4 * E0**2 - C**2)
        gamma = np.sqrt(E0**2 - (C**2 / 2))
        t1 = (E**2 - gamma**2)**2
        t2 = (0.25 * alpha**2 * C**2)
        zeta4 = t1 + t2

        # Getting epsilon(E): equation (6) in the paper
        t1 = epsR
        t2 = (A * C * a_ln) / (np.pi * zeta4 * 2 * alpha * E0) * np.log(
            (E0**2 + Eg**2 + alpha * Eg) / (E0**2 + Eg**2 - alpha * Eg))
        t3 = -1 * (A * a_atan) / (np.pi * zeta4 * E0) * (np.pi - np.arctan(
            (2 * Eg + alpha) / C) + np.arctan((-2 * Eg + alpha) / C))
        t4 = 2 * (A * E0 * Eg) / (np.pi * zeta4 * alpha) * (
            E**2 - gamma**2) * (np.pi + 2 * np.arctan(2 * (gamma**2 - Eg**2) /
                                                      (alpha * C)))
        t5 = -1 * (A * E0 * C) / (np.pi * zeta4) * (E**2 + Eg**2) / (
            E) * np.log(np.abs(E - Eg) / (E + Eg))
        t6 = 2 * (A * E0 * C) / (np.pi * zeta4) * Eg * np.log(
            (np.abs(E - Eg) *
             (E + Eg)) / (np.sqrt((E0**2 - Eg**2)**2 + Eg**2 * C**2)))
        epsR = (t1 + t2 + t3 + t4 + t5 + t6)

        # Getting imaginary part of epsilon: equation (4):
        a = 1 / E * (A * E0 * C * (E - Eg)**2) / (
            (E**2 - E0**2)**2 + (C**2) * (E**2))
        epsI = np.where(E > Eg, a, 0)

        # refractive index n:
        n = np.sqrt(0.5 * (np.sqrt(epsR**2 + epsI**2) + epsR))

        # extinction coefficient k:
        k = np.sqrt(0.5 * (np.sqrt(epsR**2 + epsI**2) - epsR))

        cls.epsilon_real = epsR
        cls.epsilon_imag = epsI
        cls.epsilon = epsR + epsI * 1j
        cls.n = n
        cls.k = k
        cls.refractive_index = n + k * 1j


def get_cplx(mat_name, wl):  # Get refractive index for a given wavelength
    lambdas = np.linspace(wl, wl)
    mat = material(mat_name, lambdas)
    return np.mean(mat.refractive_index)


def ask_for_parameters():
    while True:
        try:
            min = int(input('Minimum wavelenght: '))
            max = int(input('Maximum wavelenght: '))
            step = float(input('Stepsize: '))
            mat_name = input('Material (\'list\' for available options): ')
            break
        except:
            print('Not a valid option!')
    if mat_name == 'list':
        print(sorted(SUPPORTED_MATERIALS))
        try:
            mat_name = input('Material: ')
        except:
            print('Not a valid option!')
    return min, max, step, mat_name


def main():
    try:
        min, max, step, mat_name = ask_for_parameters()
        lambdas = np.arange(min, max + 1, step)
        mat = material(mat_name, lambdas)
        plt.figure()
        plt.plot(lambdas, mat.n, label='n')
        plt.plot(lambdas, mat.k, label='k')
        plt.legend()
        plt.show()
    except:
        print('Program closed.')


if __name__ == '__main__':  # You can load this file to quickly plot n and k.
    main()
