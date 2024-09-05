import numpy as np
import matplotlib.pyplot as plt
from cmath import *
import cmath
from matplotlib import cm
from scipy.linalg import solve_banded
import matplotlib as mpl
from scipy import special
from numpy import linalg as LA
import math
from scipy.linalg import eigh_tridiagonal
from scipy import optimize
from scipy.optimize import curve_fit
import scipy.integrate as integrate
from SE_one_particle import *

r1 = np.array([-2, 1])
r2 = np.array([-5/2, 4/3, -1/12])
r3 = np.array([-49/18, 3/2, -3/20, 1/90])
r4 = np.array([-205/72, 8/5, -1/5, 8/315,-1/560 ])
r5 = np.array([-5269/1800, 5/3, - 5/21, 5/126, - 5/1008,1/315])
r6 = np.array([-5369/1800, 12/7, -15/56, 10/189, - 1/112,2/1925, - 1/16632])
r7 = np.array([-266681/88200,7/4, - 7/24,7/108, - 7/528,7/3300, - 7/30888,1/84084])

def potencial(x, lam):
    return(0.5*x**2 + lam*x**4)

def LHO_eigen(x, N, a=0):
    Herm = special.hermite(N)
    return ((1/(np.pi**(0.25) *((2**N) * math.factorial(N))**(0.5))) *\
    Herm(x) * np.exp(-(x-a)**2.0 / 2.0))

def LHO_norm(x, N, a=0):
    lho_eigen = LHO_eigen(x,N)
    return(lho_eigen/np.sqrt(np.sum(np.abs(lho_eigen)**2)))

def overlap_integral(psi,fi,potencial):
    return(np.sum(np.conjugate(psi)*potencial*fi))

def LHO_hamilt(N, potencial,x, lambdaa):
    mat = np.zeros((N,N))
    diag = np.array([(i+0.5) for i in range(N)])
    for i in range(N):
        for j in range(N):
            if i==j:
                fja = LHO_norm(x,i)
                mat[i,i] = diag[i] + lambdaa*overlap_integral(fja,fja,potencial)
            elif j>i:
                mat[i,j] = lambdaa*overlap_integral(LHO_norm(x,i),LHO_norm(x,j),potencial)
                mat[j,i] = mat[i,j]

    return mat


def ANH_hamiltonka(N, dx, x, lam, red_r):
    V = potencial(x, lam)
    b = -1/(2.*dx**2) ### od SE 1/2
    diag_el = np.array([])
    for i in red_r:
        diag_el = np.append(diag_el, (b*i))
    d = diag_el[0] +  V
    mat = np.diag(np.ones(N)*d)
    for i in range(1,len(diag_el)):
        mat += np.diag(np.ones(N-i)*diag_el[i], k=i)
        mat += np.diag(np.ones(N-i)*diag_el[i], k=-i)
    return(mat)

#### dobro je da začetni vektor(psi0) lanczosove baze dobro aproksimira osnovno stanje. Osnovno stanje mora
#### biti normalizirano
def lanczos_metod(psi0, st_vekt_baze, hamiltonka):
    baza = np.zeros((len(psi0),st_vekt_baze))+0*1j

    zg_obdiag = np.array([])
    diagonala = np.array([])
    baza[:,0] = psi0
    d0 = (np.conjugate(psi0)@hamiltonka@psi0)

    diagonala = np.append(diagonala, d0)
    baza[:,1] = hamiltonka@psi0 - d0*psi0
    baza[:,1] = baza[:,1]/np.sqrt((np.sum(np.abs(baza[:,1])**2)))
    for i in range(2,st_vekt_baze):
        obdiag = np.conjugate(baza[:,i-2])@hamiltonka@baza[:,i-1]
        diag = np.conjugate(baza[:,i-1])@hamiltonka@baza[:,i-1]
        zg_obdiag = np.append(zg_obdiag, obdiag)
        diagonala = np.append(diagonala, diag)
        baza[:,i] = hamiltonka@baza[:,i-1]- baza[:,i-1]*diag - baza[:,i-2]*obdiag
        baza[:,i] = baza[:,i]/np.sqrt((np.sum(np.abs(baza[:,i])**2)))


    diagonala = np.append(diagonala, np.conjugate(baza[:,-1])@hamiltonka@baza[:,-1])
    zg_obdiag = np.append(zg_obdiag,np.conjugate(baza[:,-2])@hamiltonka@baza[:,-1])
    zg_obdiag = np.real(zg_obdiag)
    diagonala = np.real(diagonala)
    return baza, diagonala, zg_obdiag

def ujemanje(h,N, energije):
    #LHO_energije = np.array([0.5+i for i in range(0,N)])
    st_ujemanj = 0

    x = np.arange(-20,20,h)
    if h<0:
        return len(energije)
    else:
        dx = h
        psi0 = LHO_norm(x, 0)
        anhmat = ANH_hamiltonka(len(x), dx, x, 0.01, r2)
        baza, a,e = lanczos_metod(psi0, N, anhmat)
        #print(a, e)
        en = eigh_tridiagonal(a,e)[0]
        for i in range(len(energije)):
            razl = energije[i]-en
            if np.any(np.abs(razl) < 10**(-3)):
                st_ujemanj += 1

        neujemanja = len(en)-st_ujemanj
        #print(neujemanja)
    return(neujemanja)


def uredi_energije(energije, toleranca):
    en = np.array([energije[0]])
    k=0
    for i in range(len(energije)):
        if np.abs(energije[i]-en[k])>toleranca:
            en = np.append(en, energije[i])
            k +=1
    return en

#######################

#Nx = 4000
dx=0.05
x = np.arange(-20,20, dx)
#print(len(x))
#plt.plot(x, LHO_norm(x, 150))
#plt.show()
#dx = np.abs(x[1]-x[0])

#mat = LHO_hamilt(100,x**4,x,0.0)
#E,vektorji = LA.eigh(mat)
#print(E)


################# KALIBRACIJA NA LHO ################
lho_N = np.array([100,150,200,250,300,350,400,500,600,700,800,900,1000,1200])
lho_pravilni = np.array([9,11,14,18,21,23,25,29,32,40,45,47, 44,57])
lho_opt_dx = np.array([0.12631001980717665,0.0881966, 0.07291719223175337,0.0690983,0.06121524486553679,0.05527146,0.050628068585550516,
0.04151358852841845,0.03236068,0.03531731508884823, 0.03206583969614298, 0.03064044, 0.02685655, 0.025])


#Nar = np.linspace(100,1000,10)
#Nar = np.array([2000])


#print(ujemanje(0.017, 2000)) # 0.02 80
"""
st_ujemanj = np.array([])
minimum = np.array([])
for i in Nar:
    i= int(i)
    min = optimize.brent(ujemanje, args=(i,),brack=(0.01,0.02 ), tol=1.48e-02, maxiter=3)
    print(min)
    minimum = np.append(minimum, min)
    st_ujemanj = np.append(st_ujemanj, ujemanje(min, i))
    print(st_ujemanj)

print(st_ujemanj, minimum)
"""

def f(x, a,b):
    return a/x**b

#popt, pcov = curve_fit(f, lho_N, lho_opt_dx)
#xx = [i for i in range(100,1200,1)]
def narisi_odvisnost():
    plt.plot(xx, f(xx, *popt), label=r'$fit: 2.554/N^{0.660}$')

    plt.plot(lho_N, lho_opt_dx, 'x')
    plt.xlabel('N')
    plt.legend()
    plt.ylabel('optimalen dx')
    plt.title('LHO-Lanczos')
    plt.show()
    return None

############# 1 naloga #############
def narisi_energije():
    dx = 0.027
    N = 100
    NN = np.array([i for i in range(0,N)])
    x = np.arange(-20,20, dx)
    psi0 = LHO_norm(x, 0)
    """
    anhmat = ANH_hamiltonka(len(x), dx, x, 0, r2)
    baza, a,e = lanczos_metod(psi0, 1000, anhmat)
    en = eigh_tridiagonal(a,e)[0]
    en = uredi_energije(en, 0.5)
    LHO_energije = np.array([0.5+i for i in range(0,N)])
    mat = LHO_hamilt(150,x**4,x,0.01)
    E,vektorji = LA.eigh(mat)

    anhmat = ANH_hamiltonka(len(x), dx, x, 0.01, r2)
    baza1, a1,e1 = lanczos_metod(psi0, 1000, anhmat)
    en1 = eigh_tridiagonal(a1,e1)[0]
    en1 = uredi_energije(en1, 0.4)
    """
    mat = LHO_hamilt(150,x**4,x,0.1)
    E1,vektorji1 = LA.eigh(mat)

    anhmat = ANH_hamiltonka(len(x), dx, x, 0.1, r2)
    baza2, a2,e2 = lanczos_metod(psi0, 1000, anhmat)
    en2 = eigh_tridiagonal(a2,e2)[0]
    en2 = uredi_energije(en2, 0.2)

    mat = LHO_hamilt(150,x**4,x,1)
    E2,vektorji2 = LA.eigh(mat)

    dx = 0.026
    x = np.arange(-20,20, dx)
    psi0 = LHO_norm(x, 0)
    anhmat = ANH_hamiltonka(len(x), dx, x, 1, r2)
    baza3, a3,e3 = lanczos_metod(psi0, 1000, anhmat)
    en3 = eigh_tridiagonal(a3,e3)[0]
    en3 = uredi_energije(en3, 0.0)
    print(e3[:100])
    print(E2)


    #plt.plot(NN, LHO_energije, label=r'$\lambda=0, prave$')
    #plt.plot(NN, en[:100], label=r'$Lanczos, \lambda=0$')
    #plt.plot(NN, E[:100], label=r'$LHO \ basis, \lambda=0.01$')
    #plt.plot(NN, en1[:100], label=r'$Lanczos, \lambda=0.01$')
    plt.plot(NN, E1[:100], label=r'$LHO \ basis, \lambda=0.1$')
    plt.plot(NN, en2[:100], label=r'$Lanczos, \lambda=0.1$')
    plt.plot(NN, E2[:100], label=r'$LHO \ basis, \lambda=1$')
    plt.plot(NN, en3[:100], label=r'$Lanczos, \lambda=1$')
    plt.yscale('log')
    plt.yscale('log')
    plt.ylabel('E')
    plt.xlabel('N')
    plt.legend()
    plt.show()

#narisi_energije()


####################### konvergenca energij #########

konv_anh =np.array([3,6,7,8,11,13,17,19,23,28,30,31,32,36])
def konvergenca_lho():

    plt.plot(lho_N, konv_anh, 'x', label=r'$\lambda = 0.01$')
    plt.plot(lho_N, lho_pravilni, 'x', label=r'$\lambda =0$')
    plt.xlabel('N')
    plt.legend()
    plt.ylabel(r'$N_{kov}$')
    plt.title(r'$Lanczos$')
    plt.show()

#konvergenca_lho()
def konvergenca_anh():
    st_kon = np.array([])
    dx = 0.02
    x = np.arange(-20,20, dx)
    psi0 = LHO_norm(x, 0)
    anhmat = ANH_hamiltonka(len(x), dx, x, 0.01, r2)
    baza2, a2,e2 = lanczos_metod(psi0, 1000, anhmat)
    en2 = eigh_tridiagonal(a2,e2)[0]
    en2 = uredi_energije(en2, 0.3)
    for i in range(13,14):
        st_kon = np.append(st_kon, ujemanje(lho_opt_dx[i], lho_N[i], en2))


    return(st_kon)

#################### SEMIKLASICNA OCENA ################
def integrand(y, r, l):
    return np.sqrt(np.sqrt((r**2 + 0.25/l - y**2)/l)-0.5/l)

def razlika_ploscin(r,S_prej, l):
    Si = 4*integrate.quad(integrand, 0, r, args=(r,l))[0]
    return Si-S_prej-2*np.pi


def najdi_radije(l, Nkon):
    R = np.array([])
    S_prej = 0

    for i in range(Nkon):
        root = optimize.brentq(razlika_ploscin, 1, 200,args=(S_prej,l))
        R = np.append(R,root)
        S_prej = 4*integrate.quad(integrand, 0, root, args=(root,l))[0]

    return R

def radiji_krogov(N):
    r = np.array([np.sqrt(2)])
    for i in range(1,N):
        r = np.append(r, np.sqrt(2+r[i-1]**2))
    return r

def konv_ocena(N_zadnje_lupine, l):
    krog_max = radiji_krogov(N_zadnje_lupine)[-1]
    anh = najdi_radije(l,N_zadnje_lupine)
    maksimalen_radij_anh = 0
    for i in range(N_zadnje_lupine-1):
        if anh[i] <= krog_max and anh[i+1]> krog_max:
            maksimalen_radij_anh = anh[i]

    ocena = 4*integrate.quad(integrand, 0, maksimalen_radij_anh,
    args=(maksimalen_radij_anh,l))[0] / (np.pi * krog_max**2)

    return ocena

def narisi_oceno():
    N = np.linspace(100,1000,10)
    lambdaa = np.logspace(-2, 0, 10)
    r = np.zeros((int(len(N)), int(len(lambdaa))))
    for i in range(int(len(N))):
        for j in range(int(len(lambdaa))):
            r[i,j] = konv_ocena(int(N[i]), lambdaa[j])

    cs = plt.contourf(lambdaa, N, r, cmap='viridis')
    plt.xscale('log')
    plt.colorbar(cs, label='r')
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$N$')
    plt.show()


#narisi_oceno()



################## casovni razvoj #######################


def casovni_razvoj(t, zac_stanje, st_lastnih):

    dx = 0.02
    x = np.arange(-20,20, dx)
    psi0 = LHO_norm(x, 0)
    anhmat = ANH_hamiltonka(len(x), dx, x, 0.00, r2)
    baza, a,e = lanczos_metod(psi0, 1000, anhmat)

    E, koef = eigh_tridiagonal(a,e)

    #dxx = 0.1
    #xx = np.arange(-10,10,dxx)
    #mat = LHO_hamilt(st_lastnih,x**4,x,0.0)
    #E, koef = LA.eigh(mat)

    lastne_fje = np.zeros((len(x), st_lastnih))+ 0*1j
    for i in range(st_lastnih):
        lastne_fje[:,i] = 0
        for j in range(st_lastnih):
            lastne_fje[:,i] += koef[i,j]*baza[:,j]

    psi_t=0
    for i in range(st_lastnih):
        psi_t += lastne_fje[:,i] * np.sum(np.conjugate(lastne_fje[:,i])*zac_stanje) * np.exp(-1j * E[i]*t)
    return psi_t

def narisi_razliko(zac_st, lam):
    dt=0.04001600640256103
    t = np.arange(0,10,dt)

    dx = 0.02
    x = np.arange(-20,20, dx)
    zac_st1 = LHO_norm(x,0)
    res = resitev_2(len(x),len(t),dx,dt,x,r1,m1, zac_st1, lam )
    res1 = res[:,-1]
    res2 = res[:,0]
    res3 = res[:,int(len(t)/2)]
    #vec = casovni_razvoj(t[-1], zac_st, 100)[::10]
    #vec2 = casovni_razvoj(0, zac_st, 100)[::10]* np.sqrt(10)
    vec3 = casovni_razvoj(int(len(t)/2), zac_st, 100)

    plt.plot(x,np.real(vec3), label='lanczas')
    plt.plot(x,np.real(res3))
    #plt.plot(x,np.abs(np.real(vec)-np.real(res1)), label='t={}'.format(round(t[-1],2)))
    #plt.plot(x,np.abs(np.real(vec2)-np.real(res2)), label='t={}'.format(round(t[0],2)))
    #plt.plot(x,np.abs(np.real(vec3)-np.real(res3)), label='t={}'.format(round(t[int(len(t)/2)],2)))
    plt.xlabel('x')
    plt.ylabel(r'$|Re(\psi_{time \ evol})-Re(\psi_{propagator})|$')
    plt.title(r'$\Phi_0, \lambda=$'+'{}'.format(lam))
    #plt.yscale('log')
    plt.legend()
    plt.show()

dx = 0.02
x = np.arange(-20,20, dx)
zac_st = LHO_norm(x,0)
narisi_razliko(zac_st, 0)
