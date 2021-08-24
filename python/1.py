import numpy as np
from pyscf import gto, scf, df, dft, tools, data
from pyscf_ext import *

acfile  = '214_2019_2521_MOESM2_ESM.txt'   # HF
#acfile  = '214_2019_2521_MOESM3_ESM.txt'  # HFS
#molfile = 'H2O.xyz'
#molfile = 'O2.xyz'
molfile = 'N2.xyz'
#molfile = 'H.xyz'
basis   = 'sto3g'

###############################################################################

def charge_norm(a):
  #return 1.0 #######################
  return 0.125/(gto.gaussian_int(0, a)**3)

def read_ac(fname):
  with open(fname) as f:
    lines = f.readlines()
  basis = {'H': [[0]]}
  i=0
  while i<len(lines):
    q,nq = lines[i].split()
    q  = int(q)
    nq = int(nq)
    #print(i, '\t', q, nq)
    i+=1
    qbasis = []
    for j in range(nq):
      a,c = lines[i].split()
      a = float(a)
      c = float(c)
      #print(i, '\t', a, c)
      c = c / gto.gto_norm(0, a) * charge_norm(a)
      qbasis.append([a,c])
      i+=1
    basis[pyscf.data.elements.ELEMENTS[q]] = [[0, *qbasis]]
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
    #basis[pyscf.data.elements.ELEMENTS[q]].append( [0, [a, charge_norm(a)]] )
    basis[pyscf.data.elements.ELEMENTS[q]][0].append( [a,   charge_norm(a) / gto.gto_norm(0, a)] )
  return basis

###############################################################################

acbasis = read_ac(acfile)
add_caps(acbasis)

mol = readmol(molfile, basis=basis, ignore=True)
auxmol = df.make_auxmol(mol, acbasis)

S = mol.intor('int1e_ovlp_sph')
Hcore = mol.intor('int1e_nuc_sph') + mol.intor('int1e_kin_sph')

print(Hcore)
np.savetxt('hcore', Hcore)

overlap = {}
for q in acbasis:
    primitives = np.array(acbasis[q][0][1:])
    nprim = len(primitives)
    a = primitives[:,0]
    c = primitives[:,1]
    sums = np.zeros((nprim, nprim))
    for i in range(nprim):
        for j in range(nprim):
            sums[i,j] = a[i] + a[j]
    overlaps = gto.gaussian_int(0, sums)**3 * 8.0
    overlap[q] = c @ overlaps @ c

print(overlap)

def eri_pqi(mol, auxmol):
  pmol  = mol + auxmol
  eri3c = pmol.intor('int3c2e_sph', shls_slice=(0,mol.nbas,0,mol.nbas,mol.nbas,mol.nbas+auxmol.nbas))
  return eri3c.reshape(mol.nao_nr(), mol.nao_nr(), -1)

eri3c = eri_pqi(mol, auxmol)
#print(eri3c[0,0,0])
#print(eri3c.shape)
#print('<')
#print(eri3c)
#print('>')

Heff = np.einsum('pqi->pq', eri3c) # + Hcore

print()
print(Heff)
print()
print(acbasis['H'])
#exit(0)

#print(auxmol._basis['H'])
#print(auxmol.intor('int1e_ovlp_sph'))
#print(auxmol._basis['H'])

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
#print(mol._basis['H'])


