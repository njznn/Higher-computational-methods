import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import matplotlib.colors as colors
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy import linalg as lin
from scipy import sparse
import scipy


one_enbetadep = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/9_DMRG_II/1_betta_en_dep.txt",
skiprows=1)
def betadep():
    n = np.array([4,6,8,10,12,14,16,18,20,32])
    betasm = np.linspace(0.1,1.0,10)
    beta2 = np.linspace(1.5,9,16)
    beta = np.concatenate((betasm, beta2))

    for i in range(1,one_enbetadep[0,:].size-1):
        plt.plot(beta, -one_enbetadep[:,i], 'o--', markersize=4, label='N = {}'.format(n[i-1]))

    plt.xscale('log')
    plt.legend(loc='upper right')
    plt.ylabel(r'$E_0$')
    plt.xlabel(r'$\beta$')
    plt.grid()
    plt.show()


#betadep()
M_dep_e0 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/9_DMRG_II/1_Mdep_e0.txt",
skiprows=1)
n = np.array([14,16,18])
def Mdep():
    for i in range(1,4):
        plt.plot(M_dep_e0[:,0], -M_dep_e0[:,i], 'o--', markersize=4, label='N={}'.format(n[i-1])+',TEBD')

    plt.axhline(-M_dep_e0[0,4],label='N=14, DIAG' )
    plt.axhline(-M_dep_e0[1,4], c='orange', label='N=16, DIAG')
    plt.axhline(-M_dep_e0[2,4], c='green', label='N=18, DIAG')
    plt.ylabel(r'$E_0$')
    plt.legend()
    plt.title(r'$\beta_0 = 5$')
    plt.grid()
    plt.xlabel(r'M')

    plt.show()

#Mdep()

one_en_tebd =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/9_DMRG_II/1.energije_n_tebd.txt",
skiprows=1, usecols=(0,1))
one_en_diag =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/9_DMRG_II/1.energije_n_tebd.txt",
skiprows=1, usecols=(2), max_rows=9)

def endep():
    #plt.plot(one_en_tebd[:,0], -one_en_tebd[:,1], 'o--', markersize=4, label='TEBD')
    plt.plot(one_en_tebd[:one_en_diag.size,0], abs(one_en_diag - one_en_tebd[:one_en_diag.size,1]), 'o--', markersize=4, label='DIAG')
    plt.ylabel(r'$\frac{|E0_{TEBD}- E0_{DIAG}|}{ N}$')
    plt.xlabel(r'N')
    #plt.legend()
    plt.grid()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #plt.title(r'$\beta_0=5 \ \forall N \leq 32; \beta_0=2 \ \forall N > 32$; ')
    plt.title(r'$\beta_0=5$')
    plt.show()

#endep()
norm =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/9_DMRG_II/2_norm10.txt",
skiprows=1)
y = norm

def norm():

    t = np.linspace(0,20, 21);
    plt.plot(t, y, 'o--', markersize=4)

    plt.yscale('log')
    plt.title('dt=0.1, N=10')
    plt.xlabel('t')
    plt.ylabel(r'$|||\psi(t)||^2 - ||\psi(0)||^2|$')
    plt.grid()

    plt.show()

#norm()

#corr = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/9_DMRG_II/2_corr_80.txt",
#skiprows=1)

def correl():
    n = [i for i in range(corr[:,0].size)]
    j = [0,int(corr[:,0].size/2),corr[:,0].size/2]
    for i in range(corr[0,:].size):
        plt.plot(n, corr[:,i], 'o--', markersize=3, label="j = {}".format(j[i]) )

    plt.grid()
    plt.xlabel('i')
    plt.ylabel(r'$\langle \sigma^z_{0} \sigma^z_{j} \rangle $')
    plt.title(r'$N=80, \beta_0 = 2, M=50$')
    plt.legend()
    plt.show()

#correl()
corr = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/9_DMRG_II/2_corr_32_base_beta2.txt",
skiprows=1)

def corrtime():
    n = [i for i in range(corr[:,0].size)]
    j = [0,16,31]
    for i in range(corr[0,:].size):
        plt.plot(n, corr[:,i], 'o--', markersize=3, label="j = {}".format(j[i]) )

    plt.grid()
    plt.xlabel('t')
    plt.ylabel(r'$\langle  \ \sigma^z_{j}(t) \sigma^z_{0} \ \rangle $')
    plt.title(r'$N=32, M=50, \beta_0 = 2$')
    plt.legend()
    plt.show()

#corrtime()
ksi = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/9_DMRG_II/3_ksidep.txt",
skiprows=1)
def ksidep():

    ks = np.array([0,10**(-7),10**(-6),10**(-5),10**(-4),10**(-3),10**(-2),10**(-1)])
    n = np.array([8,10,12,14,16,18,20])
    for i in range(ksi[0,:].size):
        plt.plot(ks, ksi[:,i],'o--', markersize=3, label='N={}'.format(n[i]))


    plt.xscale('log')
    plt.grid()
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$||\max_{\{c \ \backslash \ c_{DS} \}} c||^2 / ||c_{DS}||^2$')
    plt.yscale('log')
    plt.title('M=50')
    plt.legend()
    plt.show()
#ksidep()

mag32 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/9_DMRG_II/3_3_domain_60_diffsize2.txt",
skiprows=1, usecols=[i for i in range(14)])
import matplotlib as mpl
def magnprofile():
    tsm = np.linspace(0.0,0.5,2)
    tb = np.linspace(1,11,12)

    t = np.concatenate((tsm, tb))
    #t = np.linspace(0,25,26)
    x = [i for i in range(mag32[:,0].size)]


    # Add a colorbar

    cmap = plt.get_cmap('jet',t.size*100)  # Set an empty array

    for i in range(mag32[0,:].size):
        plt.plot(x, mag32[:,i],'o--', markersize=3, color=cmap(i*100))
        print(i)

    norm = mpl.colors.Normalize(vmin=t.min(),vmax=t.max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ticks=np.linspace(t.min(),t.max(),6), label='t')

    plt.legend()
    plt.xlabel('i')
    plt.title(r'$M=50, \xi = 10^{-2}, N=60$')
    plt.ylabel(r'$\langle \sigma^z_j(t)\rangle$')

    # Display the plot
    plt.show()

magnprofile()
