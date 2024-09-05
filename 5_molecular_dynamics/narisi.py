import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import matplotlib.colors as colors
from scipy.optimize import curve_fit

tmaxdep = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/NH_tmax_T.txt",
skiprows=2)
tmaxdepmax = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/Maxwell_tmax_T.txt",
skiprows=2)

tauhoover = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/NH_tau_T.txt",
skiprows=2)
taumaxwell = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/maxwell_tau_T.txt",
skiprows=2)
dthoovers = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/NH_dt_T_same.txt",
skiprows=2)

dthoover = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/NH_dt_T.txt",
skiprows=2)
dtmaxwell = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/maxwell_dt_T.txt",
skiprows=2)

def tmaxdependance():
    indeksi = np.array([i for i in range(len(tmaxdep[:,1]))])
    """
    plt.plot(indeksi, tmaxdep[:,0], '--x', label=r'$t_{max}=1000$')
    plt.plot(indeksi, tmaxdep[:,1], '--x', label=r'$t_{max}=10000$')
    plt.plot(indeksi, tmaxdep[:,2], '--x', label=r'$t_{max}=100000$')
    plt.plot(indeksi, tmaxdep[:,3], '--x', label=r'$t_{max}=1000000$')
    """
    plt.plot(indeksi, tmaxdepmax[:,0], '--x', label=r'$t_{max}=1000$')
    plt.plot(indeksi, tmaxdepmax[:,1], '--x', label=r'$t_{max}=10000$')
    plt.plot(indeksi, tmaxdepmax[:,2], '--x', label=r'$t_{max}=100000$')
    plt.plot(indeksi, tmaxdepmax[:,3], '--x', label=r'$t_{max}=1000000$')
    ""
    plt.legend()
    plt.title(r'$Maxwell, \lambda=0, dt=0.1, t_{relax}=1000$')
    plt.xlabel('$i$')
    plt.ylabel(r'$T_{i}$')
    plt.show()

#tmaxdependance()
def dtdependance():
    indeksi = np.array([i for i in range(len(tmaxdep[:,1]))])
    """
    plt.plot(indeksi, dthoovers[:,0], '--x', label=r'$dt=0.01$')
    plt.plot(indeksi, dthoovers[:,1], '--x', label=r'$dt=0.02$')
    plt.plot(indeksi, dthoovers[:,2], '--x', label=r'$dt=0.05$')
    plt.plot(indeksi, dthoovers[:,3], '--x', label=r'$dt=0.1$')
    plt.plot(indeksi, dthoovers[:,4], '--x', label=r'$dt=0.2$')
    """
    plt.plot(indeksi, dtmaxwell[:,0], '--x', label=r'$dt=0.01$')
    plt.plot(indeksi, dtmaxwell[:,1], '--x', label=r'$dt=0.02$')
    plt.plot(indeksi, dtmaxwell[:,2], '--x', label=r'$dt=0.05$')
    plt.plot(indeksi, dtmaxwell[:,3], '--x', label=r'$dt=0.1$')
    plt.plot(indeksi, dtmaxwell[:,4], '--x', label=r'$dt=0.2$')


    plt.legend()
    plt.title(r'$maxwell, \lambda=0, t_{sample}=t_{relax}=100000$')
    plt.xlabel('$i$')
    plt.ylabel(r'$T_{i}$')
    plt.show()

#dtdependance()
def taudependance():
    indeksi = np.array([i for i in range(len(tmaxdep[:,1]))])

    """
    plt.plot(indeksi, tauhoover[:,0], '--x', label=r'$\tau=0.1$')
    plt.plot(indeksi, tauhoover[:,1], '--x', label=r'$\tau=0.5$')
    plt.plot(indeksi, tauhoover[:,2], '--x', label=r'$\tau=1.0$')
    plt.plot(indeksi, tauhoover[:,3], '--x', label=r'$\tau=1.5$')
    plt.plot(indeksi, tauhoover[:,4], '--x', label=r'$\tau=2.0$')

    """
    plt.plot(indeksi, taumaxwell[:,0], '--x', label=r'$\tau=0.2$')
    plt.plot(indeksi, taumaxwell[:,1], '--x', label=r'$\tau=0.5$')
    plt.plot(indeksi, taumaxwell[:,2], '--x', label=r'$\tau=1.0$')
    plt.plot(indeksi, taumaxwell[:,3], '--x', label=r'$\tau=1.5$')
    plt.plot(indeksi, taumaxwell[:,4], '--x', label=r'$\tau=2.0$')


    plt.legend()
    plt.title(r'$Maxwell, \lambda=0, t_{sample}=t_{relax}=100000$, dt=0.1')
    plt.xlabel('$i$')
    plt.ylabel(r'$T_{i}$')
    plt.show()

#taudependance()

#study lambda and sistem size dependance:
NHL10 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/NH_lamda_N10.txt",
skiprows=2)
NHL20 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/NH_lamda_N20.txt",
skiprows=2)
NHL30 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/NH_lamda_N30.txt",
skiprows=2)
NHL40 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/NH_lamda_N40.txt",
skiprows=2)
NHL50 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/NH_lamda_N50.txt",
skiprows=2)
NHL60 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/"\
"NH_lamda_N60.txt",skiprows=2)
NHL70 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/"\
"NH_lamda_N70.txt",skiprows=2)
NHL80 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/NH_lamda_N80.txt",
skiprows=2)
NHL100 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/"\
"NH_lamda_N100.txt",skiprows=2)

ML10 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/maxwell_lamda_N10.txt",
skiprows=2)
ML20 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/maxwell_lamda_N20.txt",
skiprows=2)

ML30 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/maxwell_lamda_N30.txt",
skiprows=2)
ML40 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/maxwell_lamda_N40.txt",
skiprows=2)
ML50 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/maxwell_lamda_N50.txt",
skiprows=2)
ML60 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/"\
"maxwell_lamda_N60.txt",skiprows=2)
ML70 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/"\
"maxwell_lamda_N70.txt",skiprows=2)
ML80 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/maxwell_lamda_N80.txt",
skiprows=2)
ML100 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/5_molecular_dynamics/"\
"maxwell_lamda_N100.txt",skiprows=2)

def lambda_dependance_T(res):
    col = plt.cm.jet(np.linspace(0,1,8))
    lamda = np.array([0,0.01,0.1, 0.2, 0.4,0.6,0.8,1])
    T = res[:,::2]

    indeksi = np.array([i for i in range(len(T[:,0]))])
    indeksi = indeksi/len(indeksi)
    for i in range(len(T[0,:])):
        plt.plot(indeksi, T[:,i],'--x' ,color = col[i], label=r'$\lambda = $'+'{}'.format(lamda[i]))
        plt.xlabel('i/N')
        plt.ylabel(r'$T_{i}$')
        plt.title('maxwell, N='+'{}'.format(len(indeksi)))
        plt.legend()
    plt.show()

#lambda_dependance_T(ML80)

def lambda_dependance_J(res):
    col = plt.cm.jet(np.linspace(0,1,8))
    lamda = np.array([0,0.01,0.1, 0.2, 0.4,0.6,0.8,1])

    T = res[:,1::2]

    indeksi = np.array([i for i in range(len(T[:,0]))])
    indeksi = indeksi/len(indeksi)
    for i in range(len(T[0,:])):
        plt.plot(indeksi, T[:,i],'--x' ,color = col[i], label=r'$\lambda = $'+'{}'.format(lamda[i]))
        plt.xlabel('i/N')
        plt.ylabel(r'$J_{i}$')
        plt.title('maxwell, N='+'{}'.format(len(indeksi)))
        plt.legend()
    plt.show()

#lambda_dependance_J(NHL10)

def N_dependance(NHL10,NHL20,NHL30,NHL40,NHL50,NHL60,NHL70,NHL80,NHL100):
    lambda_stolpec = 1 #0.1
    col = plt.cm.jet(np.linspace(0,1,9))

    plt.plot(np.array([i for i in range(len(NHL10[:,lambda_stolpec]))])\
    /(len(NHL10[:,lambda_stolpec])-1), NHL10[:,lambda_stolpec], '--x', color=col[0], label='N=10')
    plt.plot(np.array([i for i in range(len(NHL20[:,lambda_stolpec]))])\
    /(len(NHL20[:,lambda_stolpec])-1), NHL20[:,lambda_stolpec], '--x', color=col[1], label='N=20')
    plt.plot(np.array([i for i in range(len(NHL30[:,lambda_stolpec]))])\
    /(len(NHL30[:,lambda_stolpec])-1), NHL30[:,lambda_stolpec], '--x', color=col[2], label='N=30')
    plt.plot(np.array([i for i in range(len(NHL40[:,lambda_stolpec]))])\
    /(len(NHL40[:,lambda_stolpec])-1), NHL40[:,lambda_stolpec], '--x', color=col[3], label='N=40')
    plt.plot(np.array([i for i in range(len(NHL50[:,lambda_stolpec]))])\
    /(len(NHL50[:,lambda_stolpec])-1), NHL50[:,lambda_stolpec], '--x', color=col[4], label='N=50')
    plt.plot(np.array([i for i in range(len(NHL60[:,lambda_stolpec]))])\
    /(len(NHL60[:,lambda_stolpec])-1), NHL60[:,lambda_stolpec], '--x', color=col[5], label='N=60')
    plt.plot(np.array([i for i in range(len(NHL70[:,lambda_stolpec]))])\
    /(len(NHL70[:,lambda_stolpec])-1), NHL70[:,lambda_stolpec], '--x', color=col[6], label='N=70')
    plt.plot(np.array([i for i in range(len(NHL80[:,lambda_stolpec]))])\
    /(len(NHL80[:,lambda_stolpec])-1), NHL80[:,lambda_stolpec], '--x', color=col[7], label='N=80')
    plt.plot(np.array([i for i in range(len(NHL100[:,lambda_stolpec]))])\
    /(len(NHL100[:,lambda_stolpec])-1), NHL100[:,lambda_stolpec], '--x', color=col[8], label='N=100')

    plt.xlabel('i/N')
    plt.ylabel(r'$J_{i}$')
    plt.title(r'$maxwell, \lambda=1$')
    plt.legend()
    plt.show()

N_dependance(NHL10,NHL20,NHL30,NHL40,NHL50,NHL60,NHL70,NHL80,NHL100)
def func(x, a):
    return a * x

def JNN(NHL10,NHL20,NHL30,NHL40,NHL50,NHL60,NHL70,NHL80,NHL100):
    lamda = np.array([0,0.01,0.1, 0.2, 0.4,0.6,0.8,1])
    col = plt.cm.jet(np.linspace(0,1,9))
    N = np.array([10,20,30,40,50,60,70,80,100])
    vec = np.array([])
    mat = np.zeros((8, 9))
    mat[:,0] = np.sum(NHL10[:,1::2], axis=0)
    mat[:,1] = np.sum(NHL20[:,1::2], axis=0)
    mat[:,2] = np.sum(NHL30[:,1::2], axis=0)
    mat[:,3] = np.sum(NHL40[:,1::2], axis=0)
    mat[:,4] = np.sum(NHL50[:,1::2], axis=0)
    mat[:,5] = np.sum(NHL60[:,1::2], axis=0)
    mat[:,6] = np.sum(NHL70[:,1::2], axis=0)
    mat[:,7] = np.sum(NHL80[:,1::2], axis=0)
    mat[:,8] = np.sum(NHL100[:,1::2], axis=0)


    for i in range(8):
        popt, pcov = curve_fit(func, N, mat[i,:])
        vec = np.append(vec, abs(popt[0]))
        plt.plot(N, mat[i,:], '--x',color=col[i], label=r'$\lambda = $'+\
        '{}'.format(lamda[i])+ ', $\kappa$=' + '{}'.format(round(abs(popt[0]), 3)))

    plt.grid()
    plt.xlabel('N')
    plt.title('maxwell')
    plt.legend()
    plt.ylabel(r'$\langle J(N) \rangle N$')
    plt.show()

    return vec

#JNN(ML10,ML20,ML30,ML40,ML50,ML60,ML70,ML80,ML100)

lamda = np.array([0,0.01,0.1, 0.2, 0.4,0.6,0.8,1])
kappaNH = np.array([0.09237206, 0.07908758, 0.02876081, 0.01746545, 0.00998476, 0.00688947,
 0.00577459, 0.0045931])

kappaMax = np.array([0.07934786, 0.06743267 ,0.0268923,  0.01608683 ,0.00973319, 0.006999,
 0.00549507, 0.00439986])


plt.plot(lamda, kappaNH, '--x',label='NH')
plt.plot(lamda, kappaMax, '--x',label='Maxwell')
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\kappa$')
plt.legend()
#plt.show()
