"""
Script to calculate refractive indices for a range of wavelengths

-----------------------------------------
INPUT variables:

mat_name = Material name (string)
wl = Wavelength as np.linspace(min,max,points)
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

---------------------------------
Available materials for these models:

Drude-Lorentz:  Au, ITO, ITO-RTA, Ag, Al, Cu, Cr, Ni, W, Ti, Pt

-----------------
Example usage in your scripts:
import get_r_index
Au = get_r_index.material(material='Au', wl=[100,1000,1000], imSign='+')
print(Au.refractive_index)

OR:
Au = get_r_index.get_cplx(material='Au', wl=1000)      # Gets the complex refractive index for ONE wavelength
print(Au)
-------------------
You can also run this file to quickly plot n and k for

----------------------
For questions:
Nikolai Weidt: weidtn@gmail.com

"""
from math import sqrt
import numpy as np


class material():
    def __init__(self, mat_name, wl, imSign='+'):

        self.mat_name = mat_name
        self.wl = wl
        self.imSign = imSign

        if self.mat_name == 'SiO_2':

            #        # fitting coefficients (IBD)
            #        B1 = -27.424
            #        B2 = 1.3995
            #        B3 = 27.23
            #        C1 = 5643.3
            #        C2 = 10805
            #        C3 = 5407.7
            #        vali = [400 2000]
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

            # fitting coefficients  (IBD Hans)
            EgeV = [7.5554]
            AeV = [106.67]
            E0eV = [9.9221]
            CeV = [0.3209]
            eps_inf = [1.57696]
            vali = [100, 10000]

        # applied material model: Tauc-Lorentz
        #[cplx, reP, imP] = n_Tauc_Lorentz(wl, EgeV, E0eV, AeV, CeV, eps_inf, RIorEPS, imSign, vali)

        elif self.mat_name == 'SiO_2-II':

            # fitting coefficients  (IBD Taimoor)
            EgeV = [7.563]
            AeV = [77.4]
            E0eV = [9.9639]
            CeV = [0.7802]
            eps_inf = [1.72683]
            vali = [100, 10000]

        # applied material model: Tauc-Lorentz
        #[cplx, reP, imP] = n_Tauc_Lorentz(wl, EgeV, E0eV, AeV, CeV, eps_inf,                                          RIorEPS, imSign, vali)

        elif self.mat_name == 'Si_3N_4':

            # fitting coefficients (Palik, Handbook of Optical Constants of Solids, AP)
            EgeV = [4.5]
            AeV = [59.2]
            E0eV = [6.78]
            CeV = [0.49]
            eps_inf = [3.1]
            vali = [100, 10000]

            # ---- starting values
            epsR = eps_inf
            epsI = 0

            # convert wl to eV
            EeV = 1240. / wl

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
            EgeV = [3.6053]
            AeV = [111.17]
            E0eV = [4.2382]
            CeV = [2.7505]
            eps_inf = [3.32592]
            vali = [100, 10000]

        # applied material model: Tauc-Lorentz
        #[cplx, reP, imP] = n_Tauc_Lorentz(wl, EgeV, E0eV, AeV, CeV, eps_inf,                                          RIorEPS, imSign, vali)

        elif self.mat_name == 'ZrO_2-II':

            # fitting coefficients (IBD Taimoor)
            EgeV = [4.7926]
            AeV = [359.08]
            E0eV = [5.1421]
            CeV = [2.7505]
            eps_inf = [2.78023]
            vali = [100, 10000]

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

            #material model : Drude Lorentz
            self.model = 'Drude_Lorentz'
        # cplx, reP, imP = n_Drude_Lorentz(wl, w_p, amp, gam, frs, RIorEPS,                                 imSign, vali)

        elif self.mat_name == 'ITO-RTA':

            # fitting coefficients (IBD, Wetterau)
            w_p = 1.944  # plasma freq
            amp = [1, 32.248]  # oscillator strength
            gam = [0.206, 0.274]  # damping
            frs = [0, 5.136]  # resonance freq
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
            vali = [100, 10000]

            eps_inf = 0.397200
            m0 = -0.033029
            x0 = 1.012892
            k0 = 0.0241

        # applied material model: Leng-Lorentz
        ##[cplx, reP, imP] = n_Leng_Lorentz(wl, C0, Beta, Eg, Gam, Mu, eps_inf,                                          m0, x0, k0, RIorEPS, imSign, vali)

        elif self.mat_name == 'a-Si':

            # fitting coefficients (Palik, Handbook of Optical Constants of Solids, AP)
            EgeV = [1.2]
            AeV = [122]
            E0eV = [3.45]
            CeV = [2.54]
            eps_inf = [1.15]
            vali = [100, 10000]

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
            print('Error, not supported material name! Try \'list\' ')

        if self.model == 'Drude_Lorentz':
            self.n_Drude_Lorentz(wl, w_p, amp, gam, frs, imSign, vali)

#########################################################

# wl = [300,3000,10000]
# ========== dispersion model functions

# function ##[cplx, reP, imP] = n_Cauchy(lambda, C0, C1, N0, N1, N2, K0, K1, K2, RIorEPS, imSign, vali)

#           if ((min(lambda) < vali(1)) || (max(lambda) > vali(2)))
#               warning('Wavelength out ')
#               cplx = reP + 1i*imP;
#           else
#               warning('Unknown option for parameter sign');
#           end

# function ##[cplx, reP, imP] = n_Sellmeier(lambda, B1, B2, B3, C1, C2, C3, RIorEPS, imSign, vali)

#           if ((min(lambda) < vali(1)) || (max(lambda) > vali(2)))
#               warning('Wavelength out of model validity range!');
#           end

#           n_mat = sqrt(1 + B1*lambda.^2 ./ (lambda.^2 - C1) + B2*lambda.^2 ./ (lambda.^2 - C2) + B3*lambda.^2 ./ (lambda.^2 - C3));
#           k_mat = 1./(n_mat .* (B1*lambda + B2./lambda + B3./lambda.^3));

#           # defining output values
#           if strcmp(RIorEPS, 'RI')
#               reP = n_mat;
#               imP = k_mat;
#           elseif strcmp(RIorEPS, 'EPS')
#               reP = n_mat.^2 - k_mat.^2;
#               imP = 2 * n_mat .* k_mat;
#           else
#               warning('Unknown option for parameter RIorEPS');
#           end

#           if strcmp(imSign, '-')
#               cplx = reP - 1i*imP;
#           elseif strcmp(imSign, '+')
#               cplx = reP + 1i*imP;
#           else
#               warning('Unknown option for parameter sign');
#           end

    @classmethod
    def n_Drude_Lorentz(cls, wl, w_p, amp, gam, frs, imSign, vali):

        # --- for pure Drude set only amp(1) non-zero
        # --- for pure Lorentz set only amp(1) to zero
        if wl[0] < vali[0] or wl[1] > vali[1]:
            print('Wavelenght out of model validity range!')

            # convert wavelength to eV
        EeV = np.divide(1240, wl)

        epsC = 1 - np.divide(
            (amp[0] * np.square(w_p)),
            (np.square(EeV) + gam[0] * EeV * 1j))  # Drude term:

        for i in range(1, len(amp), 1):  #  Lorentz terms
            epsC += (np.divide(
                amp[i] * ((w_p)**2),
                (np.square(frs[i]) - np.square(EeV) - gam[i] * EeV * 1j)))

        cls.epsilon = epsC
        cls.epsilon_real = epsC.real
        cls.epsilon_imag = epsC.imag
        cls.n = np.sqrt(
            0.5 *
            (np.sqrt(np.square(epsC.real) + np.square(epsC.imag)) + epsC.real))
        cls.k = np.sqrt(
            0.5 *
            (np.sqrt(np.square(epsC.real) + np.square(epsC.imag)) - epsC.real))
        cls.refractive_index = cls.n + cls.k * 1j

        # function ##[cplx, reP, imP] = n_Leng_Lorentz(lambda, C0, Beta, Eg, Gam, Mu, eps_inf, m0, x0, k0, RIorEPS, imSign, vali)


#           # Source: Leng et al., JVST A, 16:3, 1654-1657 (1998)

#           if ((min(lambda) < vali(1)) || (max(lambda) > vali(2)))
#               warning('Wavelength out of model validity range!');
#           end

#           # convert lambda to eV
#           EeV = 1240./lambda;

#           epsC = eps_inf + m0.*EeV.^x0 + 1i*k0;

#           for idx = 1:1:length(C0)
#               epsC = epsC + (C0(idx)./EeV.^2) .* (exp(1i*Beta(idx)).*(Eg(idx)-EeV-1i*Gam(idx)).^Mu(idx) + exp(-1i*Beta(idx)).*(Eg(idx)+EeV+1i*Gam(idx)).^Mu(idx) - 2*real(exp(-1i*Beta(idx)).*(Eg(idx)+1i*Gam(idx)).^Mu(idx)) - 2*1i*Mu(idx).*EeV.*imag(exp(-1i*Beta(idx)).*(Eg(idx)+1i*Gam(idx)).^(Mu(idx)-1)) );
#           end

#           # defining output values
#           if strcmp(RIorEPS, 'RI')
#               reP = sqrt(0.5 * (sqrt(real(epsC).^2 + imag(epsC).^2) + real(epsC)));
#               imP = sqrt(0.5 * (sqrt(real(epsC).^2 + imag(epsC).^2) - real(epsC)));
#           elseif strcmp(RIorEPS, 'EPS')
#               reP = real(epsC);
#               imP = imag(epsC);
#           else
#               warning('Unknown option for parameter RIorEPS');
#           end

#           if strcmp(imSign, '-')
#               cplx = reP - 1i*imP;
#           elseif strcmp(imSign, '+')
#               cplx = reP + 1i*imP;
#           else
#               warning('Unknown option for parameter sign');
#           end

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

#           # defining output values
#           if strcmp(RIorEPS, 'RI')
#               reP = sqrt(0.5 * (sqrt(real(epsC).^2 + imag(epsC).^2) + real(epsC)));
#               imP = sqrt(0.5 * (sqrt(real(epsC).^2 + imag(epsC).^2) - real(epsC)));
#           elseif strcmp(RIorEPS, 'EPS')
#               reP = real(epsC);
#               imP = imag(epsC);
#           else
#               warning('Unknown option for parameter RIorEPS');
#           end

#           if strcmp(imSign, '-')
#               cplx = reP - 1i*imP;
#           elseif strcmp(imSign, '+')
#               cplx = reP + 1i*imP;
#           else
#               warning('Unknown option for parameter sign');
#           end

# function ##[cplx, reP, imP] = n_Tauc_Lorentz(lambda, EgeV, E0eV, AeV, CeV, eps_inf, RIorEPS, imSign, vali)

#           # Source: Jellison & Modine, APL 69:3, 371-373 (1996) and
#           # Erratum: Jellison & Modine, APL 69:14, 2137 (1996)

#           if ((min(lambda) < vali(1)) || (max(lambda) > vali(2)))
#               warning('Wavelength out of model validity range!');
#           end

#           # starting values
#           epsR = eps_inf;
#           epsI = 0;

#           # convert lambda to eV
#           EeV = 1240./lambda;

#           for idx = 1:1:length(AeV)

#               # ----- Absorption range (E > Egap)
#               epsI_tmp = (AeV(idx) * E0eV(idx) * CeV(idx) .* (EeV - EgeV(idx)).^2) ./ (EeV .* ((EeV.^2 - E0eV(idx)^2).^2 + (CeV(idx)^2 .* EeV.^2)));

#               # ----- Transparent range (E <= Egap)
#               zix = find(EeV <= EgeV);
#               epsI_tmp(zix) = 0;

#               # ----- Summing several oscillators
#               epsI = epsI + epsI_tmp;

#               # ----- substitutions
#               aln = (EgeV(idx)^2 - E0eV(idx)^2) .* EeV.^2 + (EgeV(idx)^2 .* CeV(idx)^2) - E0eV(idx)^2 .* (E0eV(idx)^2 + 3*EgeV(idx)^2);
#     atn = (EeV.^2 - E0eV(idx)^2) * (E0eV(idx)^2 + EgeV(idx)^2) + (EgeV(idx)^2 * CeV(idx)^2);
#     alp = sqrt(4*E0eV(idx)^2 - CeV(idx)^2);
#     gam = sqrt(E0eV(idx)^2 - 0.5*CeV(idx)^2);
#     zet = (EeV.^2 - gam^2).^2 + (0.25*alp^2 * CeV(idx)^2);

#     # ----- KKR expression from epsI (sum of several oscillators)
#     epsR = epsR + (AeV(idx) * CeV(idx) .* aln) ./ (2*pi * alp * zet * E0eV(idx)) ...
#             .* log((E0eV(idx)^2 + EgeV(idx)^2 + alp * EgeV) ./ (E0eV(idx)^2 + EgeV(idx)^2 - alp * EgeV(idx))) ...
#             - (AeV(idx) .* atn ) ./ (pi * zet * E0eV(idx)) ...
#             .* (pi - atan((2*EgeV(idx) + alp) ./ CeV(idx)) + atan((-2*EgeV(idx) + alp) ./ CeV(idx))) ...
#             + 2*AeV(idx) * E0eV(idx) * EgeV(idx) .* (EeV.^2 - gam^2) ./ (pi * zet * alp) ...
#             * (pi + 2*atan(2*(gam^2 - EgeV(idx)^2) / (alp * CeV(idx)))) ...
#             - AeV(idx) * E0eV(idx) * CeV(idx) .* (EeV.^2 + EgeV(idx)^2) ./ (pi * zet .* EeV) ...
#             .* log(abs(EeV - EgeV(idx)) ./ (EeV + EgeV(idx))) + 2*AeV(idx) * E0eV(idx) * CeV(idx) * EgeV(idx) ./ (pi * zet) ...
#             .* log(abs(EeV - EgeV(idx)) .* (EeV + EgeV(idx)) ./ sqrt((E0eV(idx)^2 - EgeV(idx)^2)^2 + EgeV(idx)^2 * CeV(idx)^2));

# end

# # defining output values
# if strcmp(RIorEPS, 'RI')
#     reP = sqrt(0.5 * (sqrt(epsR.^2 + epsI.^2) + epsR));
#     imP = sqrt(0.5 * (sqrt(epsR.^2 + epsI.^2) - epsR));
# elseif strcmp(RIorEPS, 'EPS')
#     reP = epsR;
#     imP = epsI;
# else
#     warning('Unknown option for parameter RIorEPS');
# end

# if strcmp(imSign, '-')
#     cplx = reP - 1i*imP;
# elseif strcmp(imSign, '+')
#     cplx = reP + 1i*imP;
# else
#     warning('Unknown option for parameter sign');
# end


def get_cplx(mat_name, wl):  # Get refractive index for a given wavelength
    lambdas = np.linspace(wl, wl)
    mat = material(mat_name, lambdas)
    return np.mean(mat.refractive_index)


def ask_for_parameters():
    while True:
        try:
            min = int(input('Minimum wavelenght: '))
            max = int(input('Maximum wavelenght: '))
            points = int(input('Datapoints between min and max: '))
            mat_name = input('Material (\'list\' for available options): ')
            break
        except:
            print('Not a valid option!')
    if mat_name == 'list':
        print('Au, ITO, ITO-RTA, Ag, Al, Cu, Cr, Ni, W, Ti, Pt')
        try:
            mat_name = input('Material: ')
        except:
            print('Not a valid option!')
    return min, max, points, mat_name


if __name__ == '__main__':  # You can load this file to quickly plot n vs k.
    import matplotlib.pyplot as plt
    try:
        min, max, points, mat_name = ask_for_parameters()
        lambdas = np.linspace(min, max, points)
        mat = material(mat_name, lambdas)
        plt.figure()
        plt.plot(lambdas, mat.n, lambdas, mat.k)
        plt.show()
    except:
        print('Program closed.')
