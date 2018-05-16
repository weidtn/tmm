import matplotlib.pyplot as plt
import numpy as np

lambdas, n = np.loadtxt('Popova-n.csv', delimiter=',', unpack=True)
k = np.loadtxt('Popova-k.csv')
plt.plot(lambdas,n, label='Loaded from file!')

plt.xlabel('Wavelength in nm')
plt.ylabel('n')
plt.title('Interesting Graph\nCheck it out')
plt.legend()
plt.show()
