
import numpy as np
import matplotlib.pyplot as plt
from rgpy import *
import matplotlib.cm as cm



def odvod(x,arr):
    newar = np.array([]);
    for i in range(arr.size-1):
        newar = np.append(newar, (arr[i+1]-arr[i])/(x[i+1]-x[i]))

    return newar

res = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/10_RG/heks__q=4_N12_M26_J1_new.txt",
delimiter=' ',skiprows=1)

res1 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/10_RG/potts_q4_32_J1.txt",skiprows=1)

pott = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/10_RG/potts_q4_32_J1.txt",skiprows=1)
#poten = (pott[:,1]**2-pott[:,3])/(32**2)
poten = (pott[:,1])*4/32**2
pott1 = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/10_RG/potts_16_q=2_J=05.txt",skiprows=1)
poten1 = pott1[:,0]*4/(16**2)

beta= res[:,0]
lnz = res[:,1]
lnz1 = res1[:,1]





odv = odvod(res[:,0],lnz)
odv1 = odvod(res1[:,0],odvod(res1[:,0],lnz1))
"""
fig, ax = plt.subplots()
ax_zoom = ax.inset_axes([0.1, 0.2, 0.4, 0.4])
x_zoom = beta[75:100]
y_zoom = odv[75:100]/2**16
y_zoom1 = odv1[19:26]/2**16
ax_zoom.plot(x_zoom, y_zoom)
ax_zoom.plot(res1[19:26,0], y_zoom1)
ax_zoom.plot(pott[20:27,0], poten[20:27])
"""
plt.plot(res[:-1,0], odv/3**12, label='TRG, q=5')
#plt.plot(res1[:-2,0], -odv1/np.linalg.norm(odv1), label='TRG, q=2')

plt.plot(pott[:,0], +poten, label='MC')
#plt.plot(pott1[:-1,0], -odvod(pott1[:,0],poten1)/np.linalg.norm(odvod(pott1[:,0],poten1)), label='MC')
plt.grid()
plt.ylabel(r'$C_{v}$')
plt.xlabel(r'$\beta$')

plt.legend()

plt.show()
plt.clf()

def time_cons():
    teig = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/10_RG/t_eigen_q2_B1_J05_M10",
    delimiter=' ',skiprows=0)
    tcpp = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/10_RG/t_notoptcpp_q2_B1_J05_M10",
    delimiter=' ',skiprows=0)
    tpytf = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/10_RG/time_TF_q2_B1_J0.5_M10.txt",
    delimiter=' ',skiprows=0)
    py = np.loadtxt("/home/ziga/Desktop/FMF/magisterij/visje_racunske_metode/10_RG/time_notopt_py_q2_B1_J0.5_M10.txt",
    delimiter=' ',skiprows=0)

    plt.plot(teig[:,0], teig[:,1], 'x--', label='c++, eigen')
    plt.plot(tcpp[:,0], tcpp[:,1], 'x--', label='c++, for ')
    plt.plot(tpytf[:,0], tpytf[:,1], 'x--', label='pyhton,TF')
    plt.plot(py[:5,0], py[:5,1], 'x--', label='pyhton,for')
    plt.title(r'$q=2, M=10$')
    plt.legend()
    plt.xlabel(r'$N_{iter}$')
    plt.ylabel(r'$t(s)$')
    plt.yscale('log')
    plt.show()

#time_cons()


def plot_txt_files(filenames):
  data = []
  for filename in filenames:
    data.append(np.loadtxt(filename, skiprows=1))



  plt.plot(data[0][:-1, 0], odvod(data[0][:, 0],data[0][:, 1])/3**10 ,'b',label="q=2",)
  plt.plot(data[1][:, 0], data[1][:, 1]*2/32**2, 'b--')
  plt.plot(data[2][:-1, 0], odvod(data[2][:, 0],data[2][:, 1])/3**8,'r',label="q=3")
  plt.plot(data[3][:, 0], data[3][:, 1]*2/32**2, 'r--')
  plt.plot(data[4][:-1, 0], odvod(data[4][:, 0],data[4][:, 1])/3**12,'g',label="q=4")
  plt.plot(data[5][:, 0], data[5][:, 1]*2/32**2, 'g--')
  plt.plot(data[6][:-1, 0], odvod(data[6][:, 0],data[6][:, 1])/3**12,'black',label="q=5")
  plt.plot(data[7][:, 0], data[7][:, 1]*2/32**2, '--', color='black')
  plt.plot(data[8][:-1, 0], odvod(data[8][:, 0],data[8][:, 1])/3**12,'orange',label="q=6")


  #plt.axvline(x = 1.4848, color='black', linestyle='--')
  #plt.axvline(x = 1.3169, color='black', linestyle='--')
  #plt.axvline(x = 1.6094, color='black', linestyle='--')
  #plt.axvline(x = 1.710, color='black', linestyle='--')
  #plt.axvline(x = 1.7946, color='black', linestyle='--')
  plt.xlabel(r"$\beta$")
  plt.ylabel(r"$\langle E \rangle$")
  plt.annotate('TRG', (-1.5, 2), color='red', fontsize=14,
             bbox=dict(facecolor='white', edgecolor='red'))

  plt.legend()
  plt.xlim([0,3])
  plt.grid()
  plt.text(2, -1.5, 'TRG', fontsize=12)
  plt.text(2, -1.3, ' MC', fontsize=12)
  plt.plot([1.9, 2.0], [-1.45, -1.45], color='black')
  plt.plot([1.9, 2.0], [-1.25, -1.25],"--", color='black')
  plt.title(r"$N_{iter} \geq 10, M \geq 10, J=1$")

  plt.show()

filenames = ["heks__q=2_N10_M10_J1.txt", "potts_q2_32_J1.txt","heks__q=3_N8_M9_J1.txt",
 "potts_q3_32_J1.txt", "heks__q=4_N12_M26_J1.txt","potts_q4_32_J1.txt",
 "heks__q=5_N12_M20_J1.txt","potts_q5_32_J1.txt", "heks__q=6_N12_M16_J1.txt" ]
#plot_txt_files(filenames)

def plot_matrix(filename):
  # Read the matrix from the file.
  data = np.loadtxt(filename)

  num_columns = data.shape[1] - 1
  cmap = cm.get_cmap('cool', num_columns)
  Mmax = [2,4,6,10,16,20]
  for i in range(3, data.shape[1]):
    odv = odvod(data[:,0],data[:, i])

    plt.plot(data[:-1, 0], odv/(2**(10)), color=cmap(i - 1), label=r'$M=$'+str(Mmax[i-1]))

  # Set the title and labels.
  plt.title(r"$q=6, J=0.5, N_{iter}=10$")

  #norm = plt.colors.BoundaryNorm([0, 1, 2, 3, 4], len(colo))
  #cbar = plt.colorbar(cm.ScalarMappable(norm=plt.Normalize(1, num_columns), cmap=cmap),
  #                        ticks=np.arange(1, num_columns + 1))


  #cbar.set_label(r'$N_{iter}$')

  plt.xlabel(r"$\beta$")
  plt.ylabel(r"$\langle E \rangle$")
  plt.legend()




  # Add a legend and grid.
  #plt.legend()
  plt.grid(True)

  # Show the plot.
  plt.show()

#plot_matrix("Mdep_q6_J05_Nit10.txt")
