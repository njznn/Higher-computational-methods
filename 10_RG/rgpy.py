import numpy as np
import tensorflow as tf
import scipy.linalg as la
import os
import time

os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
def q_to_decimal(q, digits):
  decimal_value = 0
  power = len(digits) - 1

  for digit in digits:
    decimal_value += digit * np.power(q, power)
    power -= 1

  return decimal_value




def twodecToBaseQ(decimalValue, q):
    digits = []
    digits.append(decimalValue // q)
    digits.append(decimalValue % q)

    return digits


def deltaf(s1, s2):
    if s1 == s2:
        return 1.0
    else:
        return 0.0

def A0_tensor(q, beta, J):
  A = np.zeros([q, q, q, q])

  for l in range(q):
    for r in range(q):
      for t in range(q):
        for d in range(q):
          A[l, r, t, d] = np.exp((beta * J / 2.0) *
                                (deltaf(t + 1, r + 1) + deltaf(l + 1, t + 1) +
                                deltaf(d + 1, l + 1) + deltaf(d + 1, r + 1)))

  ten = tf.convert_to_tensor(A)
  return ten


def ConvMat(mat: np.ndarray, F):
  if F == "F12":
    qnew = mat.shape[1]
    qold = int(np.sqrt(mat.shape[0]))
    T = np.zeros([qold, qold, qnew])
    for i in range(qnew):
      for k in range(int(np.power(qold, 2))):
        vec = twodecToBaseQ(k, qold)
        prvi = vec[0]
        drugi = vec[1]
        T[prvi, drugi, i] = mat[k, i]
  elif F == "F34":
    qnew = mat.shape[0]
    qold = int(np.sqrt(mat.shape[1]))
    T = np.zeros([qnew, qold, qold])
    for i in range(qnew):
      for k in range(int(np.power(qold, 2))):
        vec = twodecToBaseQ(k, qold)
        prvi = vec[0]
        drugi = vec[1]
        T[i, prvi, drugi] = mat[i, k]
  else:
    print("vpisi pravi parameter")
    exit(1)

  ten = tf.convert_to_tensor(T)
  return ten


def AtoFF(ten: tf.Tensor, type):
  dim = ten.shape[0]
  mat = np.zeros([int(np.power(dim, 2)), int(np.power(dim, 2))])
  row = 0
  col = 0
  if type == "sz":
    for l in range(dim):
      for d in range(dim):
        for r in range(dim):
          for t in range(dim):
            mat[q_to_decimal(dim, [l, t]), q_to_decimal(dim, [d, r])] = ten[l, r, t, d]
  elif type == "sv":
    for l in range(dim):
      for t in range(dim):
        for r in range(dim):
          for d in range(dim):
            mat[q_to_decimal(dim, [t, r]), q_to_decimal(dim, [l, d])] = ten[l, r, t, d]
  else:
    print("vpisi pravi razcep!")
    exit(1)
  return mat


def computeFF(FF: np.ndarray, maxdim):

  u, s, vh = la.svd(FF, full_matrices=False)


  mi = np.linalg.matrix_rank(FF, tol=1e-8)
  if mi%2==1:
      mi +=1;
  if mi > maxdim:
    mi = maxdim


  sqrsingvec = np.sqrt(s)
  sqrsingvec = sqrsingvec[0:mi]
  F34 = vh[0:mi, :]
  F12 = u[:, 0:mi]
  F34 = np.diag(sqrsingvec) @ F34

  F12 = F12 @ np.diag(sqrsingvec)

  F12 = -F12
  F34 = -F34

  stor = [F12, F34]

  return stor



def Anew(F3, F4, F2, F1):


  TF1 = ConvMat(F1, "F12")
  TF2 = ConvMat(F2, "F12")
  TF3 = ConvMat(F3, "F34")
  TF4 = ConvMat(F4, "F34")

  TF1 = tf.reshape(F1, [int(np.sqrt(F1.shape[0])),int(np.sqrt(F1.shape[0])),F1.shape[1]])
  TF2 = tf.reshape(F2, [int(np.sqrt(F2.shape[0])),int(np.sqrt(F2.shape[0])),F2.shape[1]])
  TF3 = tf.reshape(F3, [F3.shape[0],int(np.sqrt(F3.shape[1])),int(np.sqrt(F3.shape[1]))])
  TF4 = tf.reshape(F4, [F4.shape[0],int(np.sqrt(F4.shape[1])),int(np.sqrt(F4.shape[1]))])

  Anew = tf.einsum("lzw ,xyr,twy , zxd -> lrtd",TF3,TF1,TF4,TF2)


  return Anew

"""
def Anewcomp(lnew, rnew, tnew, dnew, F3, F4, F2, F1):

  size = int(np.sqrt(F3.shape[1]))
  sum = 0
  for x in range(size):
    for y in range(size):
      for z in range(size):
        for w in range(size):
          sum += (F3[lnew, q_to_decimal(size, [z, w])] *
           F4[tnew, q_to_decimal(size, [w, y])] * F1[q_to_decimal(size, [x, y]), rnew] *
            F2[q_to_decimal(size, [z, x]), dnew])

  return sum

def Anew(F3, F4, F2, F1):

  dim = F1.shape[1]
  An = np.zeros((dim, dim, dim, dim))
  for l in range(dim):
    for r in range(dim):
      for t in range(dim):
        for d in range(dim):
          An[l, r, t, d] = Anewcomp(l, r, t, d, F3, F4, F2, F1)
  return An

"""

def Ancontr(finaltens):
  dim = finaltens.shape[0]
  sum = tf.einsum("lltt-> ", finaltens)
  sum = tf.reduce_sum(sum)
  return sum

"""
def Ancontr(finaltens):
  dim = finaltens.shape[0]
  sum = 0
  for lr in range(dim):
    for td in range(dim):
      sum += finaltens[lr, lr, td, td]
  return sum
"""

def RGstep(Aj, maxdim):

  matsz = AtoFF(Aj, "sz")
  matsv = AtoFF(Aj, "sv")


  vecsz = computeFF(matsz, maxdim)
  vecsv = computeFF(matsv, maxdim)


  new_tensor = Anew(vecsz[1], vecsv[1], vecsv[0], vecsz[0])

  return new_tensor

def calc_norm_factor(norms):
  sum = 0
  N = len(norms)
  for i in range(N):
    sum += np.power(2, N-i) * np.log(norms[i])
  return sum

def Partitsum(numq, beta, J, steps, maxdim):
  Aj = tf.constant(A0_tensor(numq, beta, J), dtype=tf.float64)
  norms = np.zeros(steps)
  ind = 0
  while steps > 1e-4:
    Aj = RGstep(Aj, maxdim)
    norm = tf.norm(Aj)
    norms[ind] = norm
    Aj = tf.nn.l2_normalize(Aj)
    steps -= 1
    ind += 1

  res = Ancontr(Aj)
  fact = calc_norm_factor(norms)
  resall = np.array([res, fact])
  return resall

def lnZ_sweep_beta_tofile(numq, J, steps, maxdim, imedat):
  betasm = np.linspace(0.1, 1.0, 10)
  betalg = np.linspace(1.1, 3.0, 20)
  beta = np.concatenate([betasm, betalg])
  lnz = np.zeros(beta.size)
  for i in range(beta.size):
    print(((1.0*i)/(beta.size-1))*100)
    resall = Partitsum(numq, beta[i], J, steps, maxdim)
    lnz[i] = np.log(resall[0])-resall[1]

  with open(imedat, "w") as outputFile:
    outputFile.write("beta, lnZ, db=" + str(beta[1] - beta[0]) + "\n")
    for i in range(beta.size):
      outputFile.write(str(beta[i]) + " " + str(lnz[i]) + "\n")

if __name__ == "__main__":
    None

    #lnZ_sweep_beta_tofile(6,1,16,16,"lnz_16st_q=6_M16_py.txt")

    def output_numpy_arrays_to_txt_file(array1, array2, filename):
      with open(filename, "w") as f:
        for item1, item2 in zip(array1, array2):
          f.write("%s %s\n" % (item1, item2))


    def measure_time():
        koraki = np.array([2,4,6,8,10])
        casi = np.array([])
        for i in range(len(koraki)):
            if i>=5:
                casi = np.append(casi,0)
            else:
                start_time = time.time()
                print(koraki[i])
                Partitsum(2,1,1,koraki[i], 10)
                end_time = time.time()
                casi = np.append(casi,end_time - start_time )

        output_numpy_arrays_to_txt_file(koraki, casi, "time_noopt_cuda_py_q2_B1_J0.5_Nit10.txt")

    measure_time()
