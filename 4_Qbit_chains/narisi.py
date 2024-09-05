import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import matplotlib.colors as colors

N2 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/Z(B)_N=2.txt",
skiprows=2)

N4 =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/Z(B)_N=4.txt",
skiprows=2)
N8 =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/Z(B)_N=8.txt",
skiprows=2)
N12 =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/Z(B)_N=12.txt",
skiprows=2)
N12r =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/Z(B)_N=12_all_real.txt",
skiprows=2)
N16 =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/Z(B)_N=16_S2.txt",
skiprows=2)
N20 =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/Z(B)_N=20_S2.txt",
skiprows=2)

beta = np.linspace(0.01, 5, 500)

H2 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/H(B)_N=2.txt",
skiprows=2)
H4 =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/H(B)_N=4.txt",
skiprows=2)
H8 =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/H(B)_N=8.txt",
skiprows=2)
H12 =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/H(B)_N=12.txt",
skiprows=2)
H12r =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/H(B)_N=12_all_real.txt",
skiprows=2)
H16 =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/H(B)_N=16_S2.txt",
skiprows=2)
H16r =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/H(B)_N=16_S2.txt",
skiprows=2)
H20 =  np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/H(B)_N=20_S2.txt",
skiprows=2)

t = np.arange(0,5,0.01)
t1 = np.arange(0,5,0.001)

C002 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C00(t)_N=2_S5.txt",
skiprows=2)
C002r = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C00(t)_N2_S2_01_r.txt",
skiprows=2)


C0012 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C00(t)_N=12_S5.txt",
skiprows=2)

C004 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C00(t)_N=4_S5.txt",
skiprows=2)
C004r = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C00(t)_N4_S2_01_r.txt",
skiprows=2)

C008 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C00(t)_N=8_S5.txt",
skiprows=2)
C008r = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C00(t)_N8_S2_01_r.txt",
skiprows=2)

C0012r = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C00(t)_N12_S2_01_r.txt",
skiprows=2)

C0016 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C00(t)_N=16_S2_01.txt",
skiprows=2)
C0016r = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C00(t)_N16_S2_01_r.txt",
skiprows=2)

C0020 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C00(t)_N=20_S2_01.txt",
skiprows=2)
C0020r = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C00(t)_N20_S2_01_r.txt",
skiprows=2)


C224 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C22(t)_N4_S2_01_r.txt",
skiprows=2)
C228 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C22(t)_N8_S2_01_r.txt",
skiprows=2)
C2212 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C22(t)_N12_S2_01_r.txt",
skiprows=2)

C2216 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C22(t)_N16_S2_01_r.txt",
skiprows=2)

C2220 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/C22(t)_N20_S2_01_r.txt",
skiprows=2)






def makeF(bet, vector):
    return (-1/bet)*np.log(vector)

def makeH(vech, vecz):
    return(vech/vecz)



def razlike_FH():
    #plt.plot(beta,np.abs(makeF(beta,N12r[:,0]) - makeF(beta,N12r[:,1])), label='S1' )
    #plt.plot(beta,np.abs(makeF(beta,N12r[:,2]) - makeF(beta,N12r[:,1])), label='S3' )
    #plt.plot(beta,np.abs(makeF(beta,N12r[:,3]) - makeF(beta,N12r[:,1])), label='S4' )
    #plt.plot(beta,np.abs(makeF(beta,N12r[:,4]) - makeF(beta,N12r[:,1])), label='S5' )
    plt.plot(beta,np.abs(makeH(H12r[:,0],N12r[:,0]) - makeH(H12r[:,1],N12r[:,1])), label='S1' )
    plt.plot(beta,np.abs(makeH(H12r[:,2],N12r[:,2]) - makeH(H12r[:,1],N12r[:,1])), label='S3' )
    plt.plot(beta,np.abs(makeH(H12r[:,3],N12r[:,3]) - makeH(H12r[:,1],N12r[:,1])), label='S4' )
    plt.plot(beta,np.abs(makeH(H12r[:,4],N12r[:,4]) - makeH(H12r[:,1],N12r[:,1])), label='S5' )
    plt.xscale('log')
    plt.xlabel(r'$\beta$')
    plt.yscale('log')
    plt.legend()
    plt.ylabel(r'$|H(\beta)_{S2} - H(\beta)_{Si}|$')
    plt.title('N=12')

    plt.show()

#razlike_FH()

def pltF():
    #plt.plot(beta,makeF(beta, N2[:,1]), label='N=2' )
    #plt.plot(beta,makeF(beta, N4[:,1]), label='N=4' )
    #plt.plot(beta,makeF(beta, N8[:,1]), label='N=8' )
    #plt.plot(beta,makeF(beta, N12r[:,1]), label='N=12' )
    #plt.plot(beta,makeF(beta, N16), label='N=16' )
    #plt.plot(beta,makeF(beta, N20), label='N=20' )

    plt.xlabel(r'$\beta$')

    plt.ylabel(r'$F(\beta)$')
    plt.title('S2')
    plt.legend()
    plt.show()

    plt.show()

#
def pltH():

    #plt.plot(beta,makeH(H2[:,1], N2[:,1]), label='N=2' )
    #plt.plot(beta,makeH(H4[:,1], N4[:,1]), label='N=4' )
    #plt.plot(beta,makeH(H8[:,1], N8[:,1]), label='N=8' )
    plt.plot(beta,makeH(H12r[:,1], N12r[:,1]), label='N=12' )
    #plt.plot(beta,makeH(H16, N16), label='N=16' )
    #plt.plot(beta,makeH(H20, N20), label='N=20' )

    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$H(\beta)$')
    plt.title('S2')
    plt.legend()
    plt.show()

#pltH()
t1 = np.arange(0, 5,0.001 )
t = np.arange(0,5,0.01)

def narisi_C(t, t1):

    #plt.plot(t, C222, label='N=2')
    #plt.plot(t[250:],C224[250:] , label='N=4')
    #plt.plot(t[250:], C228[250:], label='N=8')
    #plt.plot(t[250:], C2212[250:], label='N=12')
    plt.plot(t[100:], C2216[100:], label='N=16 ')
    plt.plot(t[100:], C2220[100:], label='N=20')
    plt.plot(t[100:], C0016[100:], label='N=16 ')
    plt.plot(t[100:], C0020[100:], label='N=20')

    plt.xlabel('t')
    plt.legend()
    plt.ylabel(r'$C_{2 2}(t)$')
    plt.title(r'$dt=0.01, S2$')
    plt.show()

#narisi_C(t,t1)

J2 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/J(t)_N2_S2_01_r.txt",
skiprows=2)
J4 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/J(t)_N4_S2_01_r.txt",
skiprows=2)
J8 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/J(t)_N8_S2_01_r.txt",
skiprows=2)
J12 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/J(t)_N12_S2_01_r.txt",
skiprows=2)
#J12l = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/J(t)_N12_S2_01_r_long.txt",
#skiprows=2)
J16 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/4_Qbit_chains/J(t)_N16_S2_01_r.txt",
skiprows=2)










t = np.arange(0, 10,0.001 )
t1 = np.arange(0,10,0.01)



def razlika():
    plt.plot(t1, J16, label='dt=0.01, S2')
    plt.plot(t, J16S3, label='dt=0.001, S3')
    plt.xlabel('t')
    plt.ylabel(r'$\langle   \mathrm{J(t)} \cdot  \mathrm{J(0)}  \rangle$')
    plt.title('N=16')
    plt.legend()
    plt.show()

#razlika()

def Jt():

    plt.plot(t1, J2, label='N=2')
    plt.plot(t1, J4, label='N=4')
    plt.plot(t1, J8, label='N=8')
    plt.plot(t1, J12, label='N=12')
    plt.plot(t1, J16, label='N=16')




    plt.xlabel('t')
    plt.ylabel(r'$\langle   \mathrm{J(t)} \cdot  \mathrm{J(0)}  \rangle$')
    plt.title(r'$dt=0.01, S2$')
    plt.legend()
    plt.show()

Jt()

def diffusion(t, arr, N):
    integral_t = np.array([])
    for i in range(0,len(arr),1):
        integral_t = np.append(integral_t, np.trapz(arr[:i]/N, dx = np.abs(t[1]-t[0])))

    return integral_t




def plotdiff():

    plt.plot(t1, diffusion(t1,J4, 4), label='N=4')
    plt.plot(t1, diffusion(t1, J8, 8), label='N=8')
    plt.plot(t1, diffusion(t1, J12, 12), label='N=12')
    plt.plot(t1, diffusion(t1, J16, 16), label='N=16')

    plt.xlabel('t')
    plt.legend()
    plt.ylabel(r'$\int_{0}^{t} \langle   \mathrm{J(\tilde{t})} \cdot  \mathrm{J(0)}  \rangle d\tilde{t}" $')
    plt.title(r'$dt=0.01, S2$')

    plt.show()

plotdiff()


def plotdiff2():
    t_ind = np.array([200, 400, 600, 800, 999])
    res = np.zeros((5, 5))
    res[:,0] = 0
    for i in range(5):
        print(np.trapz(J4[:t_ind[i]], dx=0.01))
        res[i,1] = np.trapz(J4[:t_ind[i]], dx =0.01)
        res[i,2] = np.trapz(J8[:t_ind[i]], dx=0.01)
        res[i,3] = np.trapz(J12[:t_ind[i]], dx=0.01)
        res[i,4] = np.trapz(J16[:t_ind[i]], dx=0.01)

    NN = np.array([2,4,8,12,16])
    plt.plot(NN, res[0,:], '--x', label='T=2')
    plt.plot(NN, res[1,:], '--x', label='T=4')
    plt.plot(NN, res[2,:], '--x', label='T=6')
    plt.plot(NN, res[3,:],'--x' ,label='T=8' )
    plt.plot(NN, res[4,:],'--x' ,label='T=10' )
    plt.legend()
    plt.xlabel('dolžina verige')
    plt.ylabel(r'$\int_{0}^{t} \langle   \mathrm{J(\tilde{t})} \cdot  \mathrm{J(0)}  \rangle d\tilde{t}" $')

    plt.title(r'$dt=0.01, S2 $')

    plt.show()

plotdiff2()
