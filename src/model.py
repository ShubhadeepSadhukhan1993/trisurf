import numpy as np
import matplotlib.pyplot as plt

D=1000
beta=0.8


vmag=1.2e4

c0=1

x=np.linspace(0, 1, 100)

c=beta*vmag*np.exp(-beta*vmag*x/D)/(1-np.exp(-beta*vmag/D))/D


plt.plot(x, c)
plt.ylim(0, max(c)*1.1)
plt.xlim(0,1)
plt.show()