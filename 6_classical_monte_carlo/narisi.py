import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import matplotlib.colors as colors
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D # noqa: F401 unused import
from mpl_toolkits.mplot3d import Axes3D

pottsq216 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_16_q=2_J=05.txt",
skiprows=1)
pottsq416 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_16_q=4_J=05.txt",
skiprows=1)
pottsq516 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_16_q=5_J=05.txt",
skiprows=1)
pottsq232 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_32_q=2_J=05.txt",
skiprows=1)
pottsq432 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_32_q=4_J=05.txt",
skiprows=1)
pottsq532 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_32_q=5_J=05.txt",
skiprows=1)
pottsq264 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_64_q=2_J=05.txt",
skiprows=1)
pottsq464 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_64_q=4_J=05.txt",
skiprows=1)
pottsq564 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_64_q=5_J=05.txt",
skiprows=1)

pottsq264evenmore_relax = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_64_q=2_J=05_even_more_relax.txt",
skiprows=1)
pottsq2128 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_128_q=2_J=05.txt",
skiprows=1)
pottsq4128 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_128_q=4_J=05.txt",
skiprows=1)
pottsq5128 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_128_q=5_J=05.txt",
skiprows=1)
pottsq2128more_relax = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_128_q=2_J=05_more_relax.txt",
skiprows=1)



#pottsq2 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/potts_50_q=3_J=1.txt",
#skiprows=1)

def plot_size(ind_kol):
    beta2 = pottsq516[:,0];
    plt.plot(beta2, pottsq516[:, ind_kol]/(16**2), label='16x16');
    plt.plot(beta2, pottsq532[:, ind_kol]/(32**2), label='32x32');
    plt.plot(beta2, pottsq564[:, ind_kol]/(64**2), label='64x64');
    plt.plot(beta2, pottsq5128[:, ind_kol]/(128**2), label='128x128');
    plt.title(r'$N_{relax} =10^{8}, N_{sample}=500, J=0.5, q=5$')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\langle E \rangle / N^2$')
    #plt.ylabel(r'$\langle M \rangle / N^2$')
    plt.legend()
    plt.show()

    return None

#plot_size(1)

def mag_diff(ind_kol=2):
    beta2 = pottsq216[:,0];
    beta3 = pottsq2128more_relax[:,0];
    #plt.plot(beta2, pottsq264[:, ind_kol],'x', label=r'$64x64, N_{relax}=10^7$');
    #plt.plot(beta2, pottsq264more_relax[:, ind_kol],'x', label=r'$64x64,N_{relax}= 10^8$');
    plt.plot(beta2[::2], pottsq2128[:, ind_kol][::2],'x', label=r'$128x128, N_{relax}=10^7$');
    plt.plot(beta3, pottsq2128more_relax[:, ind_kol],'x', label=r'$128x128, N_{relax}=10^8$');
    plt.legend()
    plt.xlabel(r'$\beta$')
    #plt.ylabel(r'$\langle E \rangle / N^2$')
    plt.ylabel(r'$\langle M \rangle / N^2$')
    plt.title(r'$N_{sample}=500, J=0.5, q=2$')
    plt.show()

#mag_diff()
def suscept():
    beta2 = pottsq216[:,0];

    plt.plot(beta2,( (1/np.max((beta2)*(pottsq216[:,4]-(pottsq216[:,2]**2))))*(beta2)*(pottsq216[:,4]-(pottsq216[:,2]**2))), label=r'16x16')
    plt.plot(beta2,( (beta2/(np.max(beta2*(pottsq232[:,4]-pottsq232[:,2]**2))))*(pottsq232[:,4]-pottsq232[:,2]**2)), label=r'32x32')
    plt.plot(beta2,( (beta2/np.max((beta2)*(pottsq264[:,4]-pottsq264[:,2]**2))))*(pottsq264[:,4]-pottsq264[:,2]**2), label=r'64x64')
    plt.plot(beta2,( (beta2/np.max((beta2)*(pottsq2128[:,4]-pottsq2128[:,2]**2))))*(pottsq2128[:,4]-pottsq2128[:,2]**2), 'x', label=r'128x128')
    plt.legend()
    plt.title(r'$N_{relax}=10^8, N_{sample}=500,J=0.5, q=3 $')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\chi$')

    plt.show()

#suscept()

def cv():
    beta2 = pottsq216[:,0];
    plt.plot(beta2,( (1/np.max((beta2)*(pottsq216[:,3]-(pottsq216[:,1]**2))))*(beta2)*(pottsq216[:,3]-(pottsq216[:,1]**2))), label=r'16x16')
    plt.plot(beta2,( (beta2/(np.max(beta2*(pottsq232[:,3]-pottsq232[:,1]**2))))*(pottsq232[:,3]-pottsq232[:,1]**2)), label=r'32x32')
    plt.plot(beta2,( (beta2/np.max((beta2)*(pottsq264[:,3]-pottsq264[:,1]**2))))*(pottsq264[:,3]-pottsq264[:,1]**2), label=r'64x64')
    plt.plot(beta2,(beta2/np.max((beta2)*(pottsq2128[:,3]-pottsq2128[:,1]**2)))*abs((pottsq2128[:,3]-pottsq2128[:,1]**2)),'x', label=r'128x128')
    plt.legend()
    plt.title(r'$N_{relax}=10^7, N_{sample}=500,J=0.5, q=3 $')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$C_{v}$')
    plt.show()

#cv()

def susc_q():
    beta1 =pottsq216[:,0];
    beta2 = pottsq416[:,0];
    plt.plot(beta2, ( (beta2)*(pottsq232[:,3]-(pottsq232[:,1]**2))), label=r'q=3')
    plt.plot(beta2, ( (beta2)*(pottsq432[:,3]-(pottsq432[:,1]**2))), label=r'q=4')
    plt.plot(beta2, ( (beta2)*(pottsq532[:,3]-(pottsq532[:,1]**2))), label=r'q=5')
    plt.title(r'$N_{relax} =10^{8}, N_{sample}=500, J=0.5, lattice:32x32$')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$C_{v}$')
    #plt.ylabel(r'$\langle M \rangle / N^2$')
    plt.legend()
    plt.show()


susc_q()


def qdep(ind_kol):
    beta1 = pottsq232[:,0];
    beta2 = pottsq532[:,0];
    plt.plot(beta1, pottsq232[:, ind_kol]/(32**2), label='q=3');
    plt.plot(beta2, pottsq432[:, ind_kol]/(32**2), label='q=4');
    plt.plot(beta2, pottsq532[:, ind_kol]/(32**2), label='q=5');

    plt.title(r'$N_{relax} =10^{8}, N_{sample}=500, J=0.5, lattice:32x32$')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\langle M \rangle / N^2$')
    #plt.ylabel(r'$\langle M \rangle / N^2$')
    plt.legend()
    plt.show()

    return None


qdep(2)


# HEISENBERG:

en_size_50 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/heis_50_ENdep.txt",
skiprows=1)
en_size_500 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/heis_500_ENdep.txt",
skiprows=1)

corN100 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/heis_corN100_dep.txt",
skiprows=1)

corN100 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/heis_corN500_long_dep.txt",
skiprows=1)

def energy_of_N():
    x = en_size_500[:,0]
    plt.plot(x, en_size_500[:,1], label=r'$\beta = 0$')
    plt.plot(x, en_size_500[:,2], label=r'$\beta = 0.1$')
    plt.plot(x, en_size_500[:,3], label=r'$\beta = 1$')
    plt.plot(x, en_size_500[:,4], label=r'$\beta = 10$')
    plt.plot(x, en_size_500[:,5], label=r'$\beta = 20$')
    plt.xlabel('step')
    plt.ylabel('E')
    plt.title('N=500, J=1, h=0 ')
    plt.legend()
    plt.xscale('log')
    plt.show()

def corr():
    x = np.linspace(0,249,250)
    plt.plot(x, corN100[:,0], 'b', label=r'$\beta = 0$')
    plt.plot(-x, corN100[:,0], 'b')
    plt.plot(x, corN100[:,1], 'r',label=r'$\beta = 0.1$')
    plt.plot(-x, corN100[:,1], 'r')
    plt.plot(x, corN100[:,2], 'g',label=r'$\beta = 1$')
    plt.plot(-x, corN100[:,2], 'g')
    plt.plot(x, corN100[:,3], 'brown', label=r'$\beta = 2$')
    plt.plot(-x, corN100[:,3], 'brown')
    plt.plot(x, corN100[:,4], 'purple', label=r'$\beta = 6$')
    plt.plot(-x, corN100[:,4], 'purple')
    #plt.plot(x, corN100[:,5])
    #plt.xlim(-25,25)
    #plt.ylim(-0.1, )
    plt.grid()
    plt.xlabel('r')
    plt.ylabel('C(r)')
    plt.legend()
    plt.title(r'$N_{chain} = 500, N_{relax}=5\ 10^8, N_{sample}=5\ 10^7,h=0 $')
    plt.show()


#corr()
x = np.linspace(1,99,99)
state = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/6_classical_monte_carlo/heis_state100_b1000.txt",
skiprows=1)
fig = plt.figure()
ax = plt.axes(projection='3d')  # gca = get/create current axes
ax.quiver(x,0,0, state[:,0],state[:,1] ,state[:,2], color='k')
plt.xlim(0,100)
plt.ylim(-1,1)
ax.set_zlim(-1,1)
plt.title(r'$N_{chain} = 100, N_{relax}=1\ 10^8, N_{sample}=1\ 10^7,h=0, J=1, \beta=1000 $')

#plt.show()
