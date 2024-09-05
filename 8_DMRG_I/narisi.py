import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import matplotlib.colors as colors
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy import linalg as lin
from scipy import sparse
import scipy


err = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/8_DMRG_I/err_reconst_state.txt",
skiprows=1)
bisym = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/8_DMRG_I/entr_sym_bipart.txt",
skiprows=1)
bisymrand = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/8_DMRG_I/rand_sim_bi_entr.txt",
skiprows=1)
bientrallopc = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/8_DMRG_I/entropy_all_comp.txt",
skiprows=1)
bientrallpbc = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/8_DMRG_I/entropy_all_comp_pbc.txt",
skiprows=1)
bientrallrs = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/8_DMRG_I/entropy_all_comp_rand.txt",
skiprows=1)

ABABbi =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/8_DMRG_I/entr_ABAB_sim.txt",
skiprows=1)
AABBbi =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/8_DMRG_I/AABB_sim_entropy.txt",
skiprows=1)

n = np.array([4,6,8,10,12,14])

#reconstruction base state
"""
plt.plot(n, err[:,0], 'x--', label='OBC')
plt.plot(n, err[:,1], 'x--', label='PBC')
plt.plot(n, err[:,2], 'x--', label='RS')
plt.legend()
plt.xlabel(r'$N_{spin}$')
plt.ylabel(r'$|\psi - \psi_{MPS}|$')
plt.yscale('log')

plt.show()
"""
"""
plt.plot(n, bisym[:,0], 'x--', label='OBC')
plt.plot(n, bisym[:,1],'x--', label='PBC')
plt.xlabel('n')
plt.legend()
plt.ylabel(r'$S_{sym}$')
plt.title('Simetrična biparticija')
plt.show()
"""

def plot_bipart_entr(data):
    n = np.array([4,6,8,10,12,14])
    print(np.log(2)*2)
    for i in range(len(data[1,:])):
        plt.plot([i for i in range(1,n[i])], data[:(n[i]-1),i], 'x--', label=r'$N_{spin}$'+'={}'.format(n[i]))

    plt.legend()
    plt.title('Entropija vseh možnih biparticij sistema-RS')
    plt.ylabel('S')
    plt.xlabel('i-th bipartition')
    plt.show()


#plot_bipart_entr(bientrallrs)
"""
plt.plot(n, ABABbi[:,1], 'x--', label='PBC')
plt.plot(n, ABABbi[:,0], 'x--', label='OBC')
plt.plot(n, ABABbi[:,2], 'x--', label='RS')
plt.xlabel(r'$N_{spin}$')
plt.title('AABBAA particija')
plt.ylabel('S')
plt.legend()
plt.show()
"""
"""
plt.plot(n, bisymrand, 'x--',label='Sym. bipartition')
plt.plot(n, ABABbi[:,2],'x--' ,label='ABAB')
plt.plot(n, AABBbi[:,2], 'x--',label='AABBAA')
plt.xlabel(r'$N_{spin}$')
plt.ylabel(r'$S$')
plt.title('Entangeled entropy on random state')
plt.legend()
plt.show()
"""
