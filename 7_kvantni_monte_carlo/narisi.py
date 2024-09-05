import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import matplotlib.colors as colors
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy import linalg as lin


#beta_dep = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/1_eps_beta_M1000.txt",
#skiprows=1)
#beta_dep = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/1_eps_beta_M100.txt",
#skiprows=1)
beta_dep01 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/2_lam0.1_en_beta_M100.txt",
skiprows=1)
M_dep4 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/delez_eps_M_beta4.txt",
skiprows=1)
M_dep10 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/delez_eps_M_beta10.txt",
skiprows=1)
M_dep20 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/delez_eps_M_beta20.txt",
skiprows=1)


def plot_beta_dep():
    ep = beta_dep[:,0]
    plt.plot(ep, beta_dep[:,1], 'x--',label=r'$\lambda = 0.01$')
    plt.plot(ep, beta_dep[:,2],'x--', label=r'$\lambda = 0.1$')
    plt.plot(ep, beta_dep[:,3], 'x--',label=r'$\lambda = 0.2$')
    plt.plot(ep, beta_dep[:,4], 'x--',label=r'$\lambda = 0.4$')
    plt.plot(ep, beta_dep[:,5],'x--', label=r'$\lambda = 0.6$')
    plt.plot(ep, beta_dep[:,6],'x--', label=r'$\lambda = 0.8$')
    plt.plot(ep, beta_dep[:,7],'x--', label=r'$\lambda = 1.0$')
    plt.plot(ep, beta_dep[:,8], 'x--',label=r'$\lambda = 2.0$')
    plt.xlim(0,1.1)
    plt.xlabel(r'$\epsilon$')
    plt.ylabel(r'$delež$')
    plt.title(r'$\beta = 10,M=100, N_{relax} = 10^{6}, N_{sample} = 10^{6} $')
    plt.legend()
    plt.show()

#plot_beta_dep()


def dependance_M_epsilon(data):
    interpolated_values4 = np.array([0.40018038, 0.28345628, 0.20094653, 0.14225227,
     0.11977654, 0.10586111, 0.09518923, 0.08622298])
    interpolated_values10 = np.array([0.62596216, 0.44739474, 0.31768753, 0.22557586,
     0.18548982, 0.16129205,0.14232465, 0.13092988])

    M = np.array([50,100,200,400,600,800,1000,1200])

    interpolated_valuesb = np.array([0.06449525, 0.09508315, 0.11461371, 0.14227509,
     0.20096654, 0.28347308,0.34627741])
    interpolated_values = np.array([])
    interpolated_values100 = np.array([0.1424791, 0.20069875, 0.24721399, 0.3176002,
      0.44718806, 0.62548015, 0.75832571])

    def fja(x, a, b):
        return a*(x**b)


    beta = np.array([1,2,3,5,10,20,30])

    for i in range(2,9):
        f = interpolate.interp1d(data[:,i], data[:,0])
        interpolated_values = np.append(interpolated_values, f(0.5))


    popt1, pcov1 = curve_fit(fja, beta, interpolated_valuesb)
    popt, pcov = curve_fit(fja, beta, interpolated_values)
    popt100, pcov100 = curve_fit(fja, beta, interpolated_values100)
    plt.plot(beta, fja(beta, *popt), 'r-',
         label=r'$fit: %5.3f \ \beta^{%5.3f}$, M=1000' % tuple(popt))
    plt.plot(beta, fja(beta, *popt1), 'g-',
         label=r'$fit: %5.3f \ \beta^{%5.3f}$, M=500' % tuple(popt1))

    plt.plot(beta, fja(beta, *popt100), 'b-',
         label=r'$fit: %5.3f \ \beta^{%5.3f}$, M=100' % tuple(popt100))

    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\epsilon$')

    #plt.plot(M,interpolated_values10, 'x--',label=r'$\beta = 10$' )
    plt.plot(beta,interpolated_valuesb, 'x', color='g' )
    plt.plot(beta,interpolated_values, 'x', color='r' )
    plt.plot(beta,interpolated_values100, 'x' , color='b')
    #plt.plot(M,interpolated_values4, 'x--', label=r'$\beta = 4$')
    plt.title(r'$0.5 \ delež \  sprejema,   N_{relax} = 10^{6}, N_{sample} = 10^{6}$')

    plt.legend()
    plt.show()


en_beta = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/1_en_beta_M100.txt",
    skiprows=1)

#dependance_M_epsilon(beta_dep)
en_beta = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/1_en_beta_M100.txt",
skiprows=1)
def energy_beta(en_beta):

    betasm = np.linspace(0.3,0.7,7)
    beta = np.linspace(1,50,49)
    betall = np.concatenate((betasm, beta))
    plt.xlabel(r'$\beta$')
    plt.plot(betall, en_beta[:,0], 'x--', label=r'$\langle E \rangle$')
    plt.plot(betall, en_beta[:,1], 'x--', label=r'$\langle T \rangle$')
    plt.plot(betall, en_beta[:,2], 'x--', label=r'$\langle V \rangle$')
    plt.title(r'$M=100, N_{relax} = 10^{6}, N_{sample} = 10^{6}$')
    plt.xscale('log')
    plt.legend()
    plt.show()
en_beta01 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/2_lam0.1_en_beta_M100.txt",
skiprows=1)
en_beta02 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/2_lam0.2_en_beta_M100.txt",
skiprows=1)
en_beta05 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/2_lam0.5_en_beta_M100.txt",
skiprows=1)
en_beta1 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/2_lam1.0_en_beta_M100.txt",
skiprows=1)

def kvarenergy(en0, en1, en2,en3, en4):
    betasm = np.linspace(0.3,0.7,7)
    beta = np.linspace(1,50,49)
    betall = np.concatenate((betasm, beta))
    plt.xlabel(r'$\beta$')
    plt.plot(betall, en0[:,0], label=r'$\lambda=0$')
    plt.plot(betall, en1[:,0], label=r'$\lambda = 0.1$')
    plt.plot(betall, en2[:,0], label=r'$\lambda = 0.2$')
    plt.plot(betall, en3[:,0], label=r'$\lambda= 0.5$')
    plt.plot(betall, en4[:,0], label=r'$\lambda = 1.0$')
    plt.title(r'$M=100, N_{relax} = 10^{6}, N_{sample} = 10^{6}$')
    plt.ylim(0, 3)
    plt.ylabel(r'$\langle E  \rangle$')
    plt.legend()
    plt.show()

#kvarenergy(en_beta,en_beta01, en_beta02, en_beta05, en_beta1)

from decimal import Decimal
pi4=Decimal(np.pi**(1/4))
from numpy.polynomial.hermite import hermval
from scipy.special import hermite
from math import factorial

x = np.linspace(-5,5,1000)
def N(v):
    '''Normalization constant '''

    return 1./np.sqrt(np.sqrt(np.pi)*2**v*factorial(v))
def psi(v, x):
    """Harmonic oscillator wavefunction for level v computed on grid of points x"""

    Hr=hermite(v)

    Psix = N(v)*Hr(x)*np.exp(-0.5*x**2)

    return Psix

hist = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/lhohist_Blam1.txt",
    skiprows=1)
hist100 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/lhohist_B100.txt",
    skiprows=1)
hist1 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/lhohist_B1.txt",
    skiprows=1)

def hist():
    plt.hist(hist, bins=100,histtype='step', density=True, fill=False, label=r'$\beta=10$')

    plt.xlabel('x')
    plt.ylabel(r'$|\psi(x)|^2$')
    plt.title(r'$ N_{relax} = 10^{6}, N_{sample} = 10^{6}$')
    plt.legend()
    plt.show()
    return None



def analiticenH(N,lamb):
    qji = np.zeros((N,N))
    matrika = np.zeros((N,N))
    for i in range(N):
        matrika[i][i] = i + 1/2
        for j in range(i):
            if abs(i-j)==1:
                qji[i][j]=0.5*np.sqrt(i+j+1)
                qji[j][i]=0.5*np.sqrt(i+j+1)
    return matrika + lamb*np.linalg.matrix_power(qji,4)

Elam100 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/2_en_lam_M100.txt",
    skiprows=0, usecols=0)
Elam500 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/2_en_lam_M500.txt",
   skiprows=0, usecols=0)
Elam1000 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/7_kvantni_monte_carlo/2_en_lam_M1000.txt",
    skiprows=0, usecols=0)

def E_lambda_dep():
    lambde=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
    prejsnja = []
    for lamb in lambde:
        H=lin.eigh(analiticenH(500,lamb))[0][:100]
        prejsnja.append(np.sum(H*np.exp(-5*H)/np.sum(np.exp(-5*H))))
    plt.plot(lambde, prejsnja, 'x--', label='diagonalizacija')
    plt.plot(lambde, Elam100, 'x--', label='QMC, M=100')
    plt.plot(lambde, Elam500, 'x--',label='QMC, M=500')
    plt.plot(lambde, Elam1000, 'x--',label='QMC, M=1000')
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\langle E \rangle$')
    plt.title(r'$\beta = 100, N_{relax} = 10^{7}, N_{sample} = 10^{6}$')
    plt.legend()

    plt.show()

E_lambda_dep()
