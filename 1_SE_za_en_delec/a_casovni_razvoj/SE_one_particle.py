import numpy as np
import matplotlib.pyplot as plt
from cmath import *
import cmath
from matplotlib import cm
from scipy.linalg import solve_banded
import matplotlib as mpl

####koeficienti v krajevnem delu, r=red, prvi element:k=0, drugi element:k=1 , k je obdiagonala
r1 = np.array([-2, 1])
r2 = np.array([-5/2, 4/3, -1/12])
r3 = np.array([-49/18, 3/2, -3/20, 1/90])
r4 = np.array([-205/72, 8/5, -1/5, 8/315,-1/560 ])
r5 = np.array([-5269/1800, 5/3, - 5/21, 5/126, - 5/1008,1/315])
r6 = np.array([-5369/1800, 12/7, -15/56, 10/189, - 1/112,2/1925, - 1/16632])
r7 = np.array([-266681/88200,7/4, - 7/24,7/108, - 7/528,7/3300, - 7/30888,1/84084])
##############koeficienti v padejevi approx (casovni del), m= red, prvi element:s=1, st elementov
############## je stevilo iteracij do naslednjega koraka

m1 = np.array([-2.0])
m2 = np.array([-3.0 + 1j*1.73205, -3.0 - 1j*1.73205])
m3 = np.array([-4.64437 , -3.67781 - 1j*3.50876, -3.67781 + 1j*3.50876])
m4 = np.array([-4.20758 + 1j*5.31484, -5.79242 + 1j*1.73447, -5.79242 - 1j*1.73446, -4.20758 -1j*5.31483])
m5 = np.array([-4.64935 + 1j*7.14205, -6.70391 + 1j*3.48532, -7.29348 + 1j*0.00000, -6.70391 - 1j*3.48532, -4.64935 - 1j*7.14205])

### hermitovi polinomi
H0=1
H1= lambda x: 2*x
H2 = lambda x: 4*x**2 - 2
H3=lambda x: 8*x**3 - 12*x



def potencial(x, lam):
    return(0.5*x**2 + lam*x**4)



def H(N, dx, x, lam, red_r):
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


### z je koeficient v časovnem redu iteracij, to je končna oblika hamiltoniana
### za iteracijo A psi(t+1) = (A*) psi(t)!!

def A(N, dx, dt,x,lam, red_r, z):
    b = 1j*dt/(2*dx**2)
    V = potencial(x, lam)
    diag_el = np.array([])
    for i in red_r:
        diag_el = np.append(diag_el, (b*i)/np.conjugate(z))

    d = 1 + diag_el[0]-1j*dt*V/np.conjugate(z)
    mat = np.diag(np.ones(N)*d)
    for i in range(1,len(diag_el)):
        mat += np.diag(np.ones(N-i)*diag_el[i], k=i)
        mat += np.diag(np.ones(N-i)*diag_el[i], k=-i)

    return(mat)


def A_banded(A, red_r):
    st_diag = len(red_r)
    dim = int(len(A))
    banded = np.zeros((st_diag, dim)) + 0*1j

    for i in range(st_diag):
        banded[i,:dim-i] = np.diag(A, k=i)


    flip_mat = np.flip(banded, axis=0)
    flip_mat = np.flip(flip_mat, axis=1)
    res = np.vstack((flip_mat, banded[1:,:]))
    return res



def LHO_eigen(x, H, N, a=0):
    if N==0:
        return ((1/(np.pi**(0.25) * np.sqrt(2**N *np.math.factorial(N)))) *1 * np.exp(-0.5*(x-a)**2))
    else:
        return ((1/(np.pi**(0.25) * np.sqrt(2**N *np.math.factorial(N)))) *\
        H(x) * np.exp(-(x-a)**2 / 2))

def LHO_norm(x, H, N, a=0):
    return(LHO_eigen(x,H,N)/np.sqrt((np.sum(LHO_eigen(x,H,N)*np.conjugate(LHO_eigen(x,H,N))))))


#### koncni propagator ####

def vsota_vrste(H, psi, st_clenov, dt):
    psi = psi+0*1j
    psi_next = psi
    Hpow = H
    for i in range(1, st_clenov):

        ##racunsko neucinkovito!
        psi_next +=((-1j*dt)**i /np.math.factorial(i))*(Hpow @ psi)
        Hpow = Hpow@H
    return psi_next


def koncni_propagator(H, psi0, x,t, st_clenov):
    dx = np.abs(x[1]-x[0])
    dt = np.abs(t[1]-t[0])
    if dt<= 2*np.pi*dx**2:
        print('stabilnost ok')
    else:
        print('nestabilno, popravi razmerje dt/dx')
        exit(1)

    time_evol = np.zeros((len(x), len(t))) + 0*1j
    time_evol[:,0] = psi0
    for i in range(1, len(t)):
        time_evol[:,i] = vsota_vrste(H, time_evol[:, i-1],st_clenov,  dt)
    return(time_evol)



############### implicitna shema ###########

def padejeva(N, dx, dt,x, red_r, prejsnja_resitev, red_m, lamda):
    st_sub = len(red_r)-1
    for i in range(len(red_m)):
            mat = A(N, dx, dt,x, lamda,red_r, red_m[i])
            Hconj = np.conjugate(mat)
            matrika_band = A_banded(mat, red_r)
            res = solve_banded((st_sub,st_sub),matrika_band, Hconj@ prejsnja_resitev)
            prejsnja=res
    return(res)

def resitev_1(zacetni_pogoj, N, n,matrika, red_r):
    st_sub = len(red_r)-1
    res = np.zeros((N, n), dtype=complex)
    res[:,0] = zacetni_pogoj
    H_banded = A_banded(matrika, red_r)
    mat = matrika
    Hconj = np.conjugate(mat)
    for i in range(1, n):
            res[:, i] = solve_banded((st_sub,st_sub),H_banded, Hconj@(res[:, i-1]))
    return res

def resitev_2(N, n, dx, dt, x, red_r, red_m, zacetni_pogoj, lamda):
    st_sub = len(red_r)-1
    res = np.zeros((N, n), dtype=complex)
    res[:,0] = zacetni_pogoj
    for i in range(1, n):
            res[:, i] = padejeva(N, dx, dt,x, red_r, res[:, i-1],red_m, lamda)
    return res




#####################################
"""
Nx = 200
x = np.linspace(-10,10,Nx)
dx = np.abs(x[1]-x[0])
#dt=t[1]-t[0]
dt=0.004001600640256103
t = np.arange(0,20,dt)
Nt = len(t)
#print(Nt)

dt=t[1]-t[0]
#print(dt, 2*np.pi*dx**2)

zac_pog = LHO_eigen(x, H1, 1)

print(np.sum(LHO_norm(x,H0,0)*np.conjugate(LHO_norm(x,H1,1))))
"""
############### grafi ############

def pdf(res):
    pdft = np.array([])
    for i in range(len(res[0,:])):
        vsota = 0
        for x in range(len(res[:,0])):
            vsota += np.abs(res[x,i])**2
        pdft = np.append(pdft, vsota)

    return pdft/np.sum(np.abs(res[:,0])**2)


#lambdaa = 0.1
#Hamilt = H(Nx, dx, x, lambdaa, r2)
#print(pdf(koncni_propagator(Hamilt, zac_pog, x, t,10 )))
#matrika_sistema  = A(Nx, dx, dt,x,lambdaa, r1, m1[0])
#res1 = resitev_1(zac_pog,Nx, Nt,matrika_sistema,  r1 )
#res2= resitev_2(Nx, Nt,dx, dt, x, r2, m1, zac_pog, lambdaa)
#res3 = koncni_propagator(Hamilt, zac_pog, x, t,10 )
#print(res3)
#Hamilt = H(Nx, dx, x, lambdaa, r2)
#res4 = koncni_propagator(Hamilt, zac_pog, x, t,10 )


def wave_funciton(x):
    plt.plot(x, LHO_eigen(x, H0,0), label=r'$\phi_0$')
    plt.plot(x, np.abs(LHO_eigen(x, H1, 1)), label=r'$\phi_1$')
    plt.plot(x, np.abs(LHO_eigen(x, H2, 2)), label=r'$\phi_2$')
    plt.yscale('log')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel(r'$|\psi(x)|$')
    plt.title(r'$t=0, \lambda = 0$')
    plt.show()

#wave_funciton(x)
def plot_pdf(t, res1, res2, res3, res4):
    plt.plot(t, pdf(res1), label='implicit, m=1,r=1')
    plt.plot(t, pdf(res2),  label='implicit, m=1,r=2')
    plt.plot(t, pdf(res3),  label='jump,r=1')
    plt.plot(t, pdf(res4),  label='jump,r=2')
    plt.legend()
    plt.xlabel('t')
    plt.ylabel(r'$|\phi_0(t)|^2$')
    plt.title(r'$\lambda = 0$')
    plt.show()

#plot_pdf(t,res1, res2, res3, res4 )

def time_evol(t,x, res1, res2,res3, res4):
    fig, ax = plt.subplots(2,2)
    norm = mpl.colors.Normalize(vmin=t[0],vmax=t[-1])
    for i in range(0,len(t),100):
        cmap = plt.get_cmap('viridis',len(t))
        ax[0][0].plot(x, np.real(res1[:,i]),c=cmap(i))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    plt.colorbar(sm, ax=ax[0][0], label='t')
    ax[0][0].set_xlabel('x')
    ax[0][0].set_title(r'$\phi_0$')
    ax[0][0].set_ylabel(r'$Re(\psi)$')

    for i in range(0,len(t),100):
        cmap = plt.get_cmap('viridis',len(t))
        ax[0][1].plot(x, np.real(res2[:,i]),c=cmap(i))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    plt.colorbar(sm, ax=ax[0][1], label='t')
    ax[0][1].set_xlabel('x')
    ax[0][1].set_title(r'$\phi_1$')
    ax[0][1].set_ylabel(r'$Re(\psi)$')

    for i in range(0,len(t),100):
        cmap = plt.get_cmap('viridis',len(t))
        ax[1][0].plot(x, np.real(res3[:,i]),c=cmap(i))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    plt.colorbar(sm, ax=ax[1][0], label='t')
    ax[1][0].set_xlabel('x')
    ax[1][0].set_title(r'$\phi_2$')
    ax[1][0].set_ylabel(r'$Re(\psi)$')

    for i in range(0,len(t),100):
        cmap = plt.get_cmap('viridis',len(t))
        ax[1][1].plot(x, np.real(res4[:,i]),c=cmap(i))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    plt.colorbar(sm, ax=ax[1][1], label='t')
    ax[1][1].set_xlabel('x')
    ax[1][1].set_title(r'$\phi_3$')
    ax[1][1].set_ylabel(r'$Re(\psi)$')

    plt.suptitle(r'$\lambda=0.5$')
    plt.tight_layout(pad=0.01)
    plt.show()


"""
zac_pog = LHO_eigen(x, H0, 0)
res1 = resitev_2(Nx,Nt,dx,dt,x,r3,m2, zac_pog, lambdaa )
zac_pog = LHO_eigen(x, H1, 1)
res2 = resitev_2(Nx,Nt,dx,dt,x,r3,m2, zac_pog, lambdaa )
zac_pog = LHO_eigen(x, H2, 2)
res3 = resitev_2(Nx,Nt,dx,dt,x,r3,m2, zac_pog, lambdaa )
zac_pog = LHO_eigen(x, H3, 3)
res4 = resitev_2(Nx,Nt,dx,dt,x,r3,m2, zac_pog, lambdaa  )
"""

#time_evol(t, x, res3/np.sum(np.abs(res3[:,0])**2), res3/np.sum(np.abs(res3[:,0])**2), \
#res3/np.sum(np.abs(res3[:,0])**2),res3/np.sum(np.abs(res3[:,0])**2) )

#zac_pog = LHO_eigen(x, H1, 1)
def odstopanja(zac_pog_,t_, l_, Nx_,Nt_ ,dx_, dt_, x_):
    referencna =resitev_2(Nx_,Nt_,dx_,dt_,x_,r4,m4, zac_pog_, l_)
    res1 = resitev_2(Nx_,Nt_,dx_,dt_,x_,r1,m1, zac_pog_, l_ )
    Hamilt = H(Nx_, dx_, x_, l_, r1)
    res2 = koncni_propagator(Hamilt, zac_pog_, x_, t_,10 )
    Hamilt = H(Nx_, dx_, x_, l_, r2)
    res3 = koncni_propagator(Hamilt, zac_pog_, x_, t_,10 )

    plt.plot(x, np.abs(np.abs(referencna[:,-1])-np.abs(res1[:,-1])), label='implicit, r=1,m=1')
    plt.plot(x, np.abs(np.abs(referencna[:,-1])-np.abs(res2[:,-1])), label='jump, r=1')
    plt.plot(x, np.abs(np.abs(referencna[:,-1])-np.abs(res3[:,-1])), label='jump, r=2')
    plt.yscale('log')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel(r'$||\psi_{ref}| - |\psi_i||$')
    plt.title(r'$\phi_1, t=10, \lambda=0.1$')
    plt.show()

#odstopanja(zac_pog, t, 0.1, Nx,Nt, dx, dt, x)


############ 2 naloga ################
"""
lambdaa = 0.1
zac_pog = LHO_eigen(x, H0, 0,1.2)
res1 = resitev_2(Nx,Nt,dx,dt,x,r3,m2, zac_pog, 0.01 )
zac_pog = LHO_eigen(x, H0, 0,1.2)
res2 = resitev_2(Nx,Nt,dx,dt,x,r3,m2, zac_pog, 0.03 )
zac_pog = LHO_eigen(x, H0, 0,1.2)
res3 = resitev_2(Nx,Nt,dx,dt,x,r3,m2, zac_pog, 0.06 )
zac_pog = LHO_eigen(x, H0, 0,1.2)
res4 = resitev_2(Nx,Nt,dx,dt,x,r3,m2, zac_pog, 0.09 )
"""
def coherent_evol(t, x,res1, res2, res3, res4):
    fig, ax = plt.subplots(2,2)
    norm = mpl.colors.Normalize(vmin=t[0],vmax=t[-1])
    for i in range(0,len(t),1000):
        cmap = plt.get_cmap('viridis',len(t))
        ax[0][0].plot(x, np.real(res1[:,i]),c=cmap(i))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    plt.colorbar(sm, ax=ax[0][0], label='t')
    ax[0][0].set_xlabel('x')
    ax[0][0].set_title(r'$\lambda=0.01$')
    ax[0][0].set_ylabel(r'$Re(\psi)$')

    for i in range(0,len(t),1000):
        cmap = plt.get_cmap('viridis',len(t))
        ax[0][1].plot(x, np.real(res2[:,i]),c=cmap(i))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    plt.colorbar(sm, ax=ax[0][1], label='t')
    ax[0][1].set_xlabel('x')
    ax[0][1].set_title(r'$\lambda=0.03$')
    ax[0][1].set_ylabel(r'$Re(\psi)$')

    for i in range(0,len(t),1000):
        cmap = plt.get_cmap('viridis',len(t))
        ax[1][0].plot(x, np.real(res3[:,i]),c=cmap(i))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    plt.colorbar(sm, ax=ax[1][0], label='t')
    ax[1][0].set_xlabel('x')
    ax[1][0].set_title(r'$\lambda=0.06$')
    ax[1][0].set_ylabel(r'$Re(\psi)$')

    for i in range(0,len(t),1000):
        cmap = plt.get_cmap('viridis',len(t))
        ax[1][1].plot(x, np.real(res4[:,i]),c=cmap(i))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    plt.colorbar(sm, ax=ax[1][1], label='t')
    ax[1][1].set_xlabel('x')
    ax[1][1].set_title(r'$\lambda=0.09$')
    ax[1][1].set_ylabel(r'$Re(\psi)$')

    plt.suptitle(r'$a=1.2$')
    plt.tight_layout(pad=0.01)
    plt.show()

#coherent_evol(t, x, res1/np.sum(np.abs(res1[:,0])**2), res2/np.sum(np.abs(res2[:,0])**2), \
#res3/np.sum(np.abs(res3[:,0])**2),res4/np.sum(np.abs(res4[:,0])**2) )
