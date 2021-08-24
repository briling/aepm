import numpy as np
from pyscf import gto, scf, df, dft, tools, data
from pyscf_ext import *

acfile  = '214_2019_2521_MOESM2_ESM.txt'   # HF
#acfile  = '214_2019_2521_MOESM3_ESM.txt'  # HFS
#molfile = 'H2O.xyz'
#molfile = 'O2.xyz'
molfile = 'HO.xyz'
#molfile = 'N2.xyz'
#molfile = 'H.xyz'
basis   = 'sto3g'

###############################################################################

def renormalize(a):
  # 1/norm1 = \int \exp(-a*r^2) d^3 r       => norm1 = (a/pi)^(3/2)
  # 1/norm2^2 = \int (\exp(-a*r^2))^2 d^3 r => norm2 = (2.0*a/pi)^(3/4)
  # coefficient = norm1 / norm2 = (0.5*a/pi)^(3/4)
  x = np.sqrt(np.sqrt(0.5*a/np.pi))
  return x*x*x

def read_ac(fname):
  with open(fname) as f:
    lines = f.readlines()
  basis = {'H': []}
  il=0
  while il<len(lines):
    q,ng = map(int,lines[il].split())
    il+=1
    qbasis = []
    for ig in range(ng):
      a,c = map(float,lines[il].split())
      qbasis.append([0,[a, c*renormalize(a)]])
      il+=1
    basis[pyscf.data.elements.ELEMENTS[q]] = qbasis
  return basis

def add_caps(basis):
  caps_array = np.zeros(103)
  caps_array  [  1 :   2 +1] = 1.0 /  3.0,
  caps_array  [  3 :   4 +1] = 1.0 / 16.0,
  caps_array  [  5 :  10 +1] = 1.0 /  3.0,
  caps_array  [ 11 :  12 +1] = 1.0 / 32.0,
  caps_array  [ 13 :  18 +1] = 1.0 /  8.0,
  caps_array  [ 19 :  20 +1] = 1.0 / 32.0,
  caps_array  [ 21 :  30 +1] = 1.0 /  6.0,
  caps_array  [ 31 :  36 +1] = 1.0 / 12.0,
  caps_array  [ 37 :  38 +1] = 1.0 / 32.0,
  caps_array  [ 39 :  48 +1] = 1.0 /  8.0,
  caps_array  [ 49 :  54 +1] = 1.0 / 12.0,
  caps_array  [ 55 :  70 +1] = 1.0 / 32.0,
  caps_array  [ 71 :  86 +1] = 1.0 / 12.0,
  caps_array  [ 87 : 102 +1] = 1.0 / 32.0
  for q in range(1,103):
    a = caps_array[q]
    basis[pyscf.data.elements.ELEMENTS[q]].append( [0, [a, renormalize(a) ]] )
  return basis

def Heff_LB2020(mol, acbasis):
  auxmol = df.make_auxmol(mol, acbasis)
  pmol  = mol + auxmol
  eri3c = pmol.intor('int3c2e_sph', shls_slice=(0,mol.nbas,0,mol.nbas,mol.nbas,mol.nbas+auxmol.nbas))
  eri3c = eri3c.reshape(mol.nao_nr(), mol.nao_nr(), -1)
  iao = 0
  for iat in range(auxmol.natm):
    q = auxmol._atom[iat][0]
    for prim in auxmol._basis[q]:
      eri3c[:,:,iao] *= prim[1][1]
      iao+=1
  return np.einsum('pqi->pq', eri3c)


###############################################################################

acbasis = read_ac(acfile)
add_caps(acbasis)

mol = readmol(molfile, basis=basis, ignore=True)

S = mol.intor('int1e_ovlp_sph')
Hcore = mol.intor('int1e_nuc_sph') + mol.intor('int1e_kin_sph')

Heff = Hcore + Heff_LB2020(mol, acbasis)

print()
print(Heff)
print()
# print(acbasis['H'])
# #exit(0)

#a = 1
#print ( a * np.sqrt(a) / (np.pi * np.sqrt(np.pi))  )
#print ( 0.125/(gto.gaussian_int(0, a)**3)  )

###############################################################################
###############################################################################
###############################################################################
###############################################################################





#np.savetxt('nuc.dat', mol.intor('int1e_kin_sph'))
#np.savetxt('core.dat', Hcore)
#print(Hcore[0,0])
#print(mol.intor('int1e_nuc_sph')[0,0])
#
#for i in mol._basis['O']:
#  for j in i:
#    print(type(j))
#  print()
#

#a = 1.8
#print(np.power(0.5*a*np.pi, 0.25))
#print(np.power(0.5*a*np.pi, 3))
