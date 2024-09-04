import numpy as np
import matplotlib.pyplot as plt
from cmath import *
from matplotlib import cm
from scipy.stats import moment

############### 1 paket ##############
#fizikalne kostnante-globale spremenljivke:
w=1
lam = 0
k = w**2
alp = k**(1/4)

#mreža:
N = 300
n = 10000 #velikost mreže v času
T = 2*np.pi/w
x = np.linspace(-40,40,N)
t = np.linspace(0, 10*T, n)
dt = t[1]-t[0]
dx = x[1]-x[0]
print(dx**2)
print(dt)


def koherentno(x, ti):
    ksi = alp*x
    ksi2 = alp*lam
    res = np.sqrt(alp/np.sqrt(np.pi)) * np.exp(-0.5*(ksi-ksi2*np.cos(w*ti))**2 -
    1j*((w*ti/2) + ksi*ksi2*np.sin(w*ti) - 0.25*ksi2**2 * np.sin(2*w*ti)))
    return res

def zacetni_pogoj1(x):
    return np.sqrt(alp/np.sqrt(np.pi)) * np.exp(-0.5*alp**2*(x-lam)**2)

########## 2. paket ################
sigma0 = 1/20
k0 = 50*np.pi
lamb =0.25
N2 = 600
x2 = np.linspace(-0.5, 1.5, N2)
dx2 = x2[1] - x2[0]
dt2 = 2*dx2**2
t2 = np.arange(0,0.0032 + dt2, dt2 )

def zacetni_gauss(x):
    return (2*np.pi*sigma0**2)**(-1/4)* np.exp(1j*k0*(x-lamb))* np.exp(-(x-lamb)**2/(2*sigma0)**2)

def analiticna_gauss(x, ti):
    predfaktor = ((2*np.pi*sigma0**2)**(-1/4))/np.sqrt(1 + 1j*ti/(2*sigma0**2))
    eks = np.exp((-(x-lamb)**2/(2*sigma0)**2 + 1j*k0*(x-lamb) - 1j*k0**2*ti/2)/(1 + 1j*ti/(2*sigma0**2)))
    return predfaktor*eks



#########dolocitev casa, ko je tezisce x=0.75#######
def ocena_casa():
    for i in t2:
        abs = np.abs(analiticna_gauss(x2,i)**2)
        if np.abs(x2[np.where(abs == np.max(abs))] - 0.7500) < 1e-3:
            #print(x2[np.where(abs == np.max(abs))])
            print(i)
            break




if __name__ == "__main__":

    def plot_surface():
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        X,Y = np.meshgrid(t[::100],x, sparse=False)
        matrika = np.zeros((len(x), 100), dtype=complex )
        for j, i in enumerate(t[::100]):
            matrika[:,j] = koherentno(x, i)
            print(matrika[:,j])


        surf = ax.contour3D(X, Y,np.abs(matrika)**2, 50, cmap='binary',
                               linewidth=0.01, antialiased=False)
        fig.tight_layout()
        plt.show()

    plot_surface()

    def ob_periodah():
        for i in range(0,1000, 100):

            plt.plot(x, np.abs(koherentno(x,0))**2 - np.abs(koherentno(x, t[i]))**2, label='{}'.format(i/1000))
        plt.legend()
        plt.show()
    #ob_periodah()
    def plot_surface_gauss():
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        X,Y = np.meshgrid(t2,x2, sparse=False)
        matrika = np.zeros((N2, len(t2)), dtype=complex )
        for j, i in enumerate(t2):
            matrika[:,j] = analiticna_gauss(x2, i)


        surf = ax.contour3D(X, Y,np.abs(matrika)**2,1000, cmap='viridis',
                               linewidth=0, antialiased=False)
        fig.tight_layout()
        plt.show()
    #plot_surface_gauss()
    t3 = np.linspace(0.0001,0.003, 100)
    for i in t3:
        plt.plot(x2, np.abs(analiticna_gauss(x2, i))**2)
    #plt.show()
