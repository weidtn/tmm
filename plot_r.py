import matplotlib.pyplot as plt
import get_r_index
import numpy as np
from numpy import linspace

# material = input('Material name: ')
wavelength = [300, 3000, 10000]
# RIorEPS = input('RI or EPS: ')
# imSign = input('n+jk or n-jk: ')

material = 'Al'
RIorEPS = 'RI'
imSign = '+'

cplx, reP, imP = get_r_index.get_r_index(material, wavelength, RIorEPS, imSign)
lambdas = np.asarray(linspace(wavelength[0], wavelength[1], wavelength[2]))
plt.figure()
plt.plot(lambdas, cplx.real, lambdas, cplx.imag)
plt.xlabel('wavelength')
plt.ylabel('n and k')
plt.show()
