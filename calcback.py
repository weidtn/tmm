import numpy as np
import get_r_index
import matplotlib.pyplot as plt
import pandas as pd

def init(df):

    df2 = df 
    print(df2)

# LAMBDAS = np.linspace(300,2400)
LAMBDAS = np.arange(300,2700+1) 
Psi = np.array(128)
Delta = np.array(43)
df = pd.DataFrame({'lambda':LAMBDAS, 'n_S':np.round(get_r_index.material('a-Si',LAMBDAS).refractive_index,4), 'Psi':Psi, 'Delta':Delta})
df.set_index('lambda', inplace=True)
phi_i = np.deg2rad(70)  # incident angle
d_L = 100 # layer thickness
step = 0.001 # step size
n_air = 1 # Refractive index of air
n_S = 3.6499 # refractive index of substrate
rho_giv = np.tan(Psi) * np.exp(1j*Delta)
start = 0.1 # initial refractive index
n_L = start
data = np.empty((3,8000)) # array [[wl], [diff], [n_L]]
lambda_vac = LAMBDAS[1000]

# for lambda_vac, row in df.iterrows():
#     print(lambda_vac)
for i in range (1, 8000):

    # Snells Law:
    phi_L = np.arcsin(np.deg2rad((np.sin(phi_i)*n_air)/n_L))
    # phi_S = np.arcsin((np.sin(phi_L)*n_L)/n_S)
    phi_S = np.arcsin(np.deg2rad((np.sin(phi_i)*n_air)/n_S))

    # Fresnel equations:
    # air/layer:
    rs_al = (n_air*np.cos(phi_i) - n_L * np.cos(phi_L)) / n_air * np.cos(phi_i + n_L * np.cos(phi_L))
    rp_al = (n_L * np.cos(phi_i) - n_air * np.cos(phi_L)) / n_L * np.cos(phi_i) + n_air * np.cos(phi_L)

    # layer/substrate:
    rs_ls = n_L * np.cos(phi_L) - n_S * np.cos(phi_S) / n_L * np.cos(phi_L) + n_S * np.cos(phi_S)
    rp_ls = n_S * np.cos(phi_L) - n_L * np.cos(phi_S) / n_S * np.cos(phi_L) + n_L * np.cos(phi_S)


    beta = (2*np.pi/lambda_vac) * d_L * n_L * np.cos(phi_L)

    rp_L = (rp_al + rp_ls * np.exp(1j * 2 * beta)) / (1 + rp_al * rp_ls * np.exp(1j * 2 * beta))
    rs_L = (rs_al + rs_ls * np.exp(1j*2*beta)) / (1+ rs_al * rs_ls * np.exp(1j * 2 * beta))

    rho_L = rp_L / rs_L
    data[0][i] = lambda_vac
    data[1][i] = abs(rho_L - rho_giv) 
    data[2][i] = n_L
    n_L = np.round(n_L + step, 4)
# print(data)
# print(np.argmin(diff))
# print(start + (np.argmin(diff)*step))
# print(diff[np.argmin(diff)])
# plt.plot(data[1])
# plt.show()
