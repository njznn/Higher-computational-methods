
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import matplotlib.colors as colors

data = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/3_T-S_and_simplectic_integration/razlike.txt",
skiprows=1)
print(data[:,0][-1])

def razlike():
    plt.plot(data[:,0][::1000], data[:,1][::1000])
    plt.plot(data[:,0][::1000], data[:,2][::1000])
    plt.plot(data[:,0][::1000], data[:,3][::1000])
    plt.plot(data[:,0][::1000], data[:,4][::1000])
    plt.yscale('log')
    plt.show()

razlike()
