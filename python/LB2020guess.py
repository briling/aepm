import numpy
import pyscf.data, pyscf.df

class LB2020guess:

  acfile_default  = '../data/parameters_HF.dat'

  def __init__(self, fname=None):
    self.get_basis(fname)

  def renormalize(self, a):
    # 1/norm1 = \int \exp(-a*r^2) d^3 r       => norm1 = (a/pi)^(3/2)
    # 1/norm2^2 = \int (\exp(-a*r^2))^2 d^3 r => norm2 = (2.0*a/pi)^(3/4)
    # coefficient = norm1 / norm2 = (0.5*a/pi)^(3/4)
    x = numpy.sqrt(numpy.sqrt(0.5*a/numpy.pi))
    return x*x*x

  def read_ac(self, fname):
    if fname==None:
      fname = self.acfile_default
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
        qbasis.append([0,[a, c*self.renormalize(a)]])
        il+=1
      basis[pyscf.data.elements.ELEMENTS[q]] = qbasis
    return basis

  def add_caps(self, basis):
    caps_array = numpy.zeros(103)
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
      basis[pyscf.data.elements.ELEMENTS[q]].append( [0, [a, self.renormalize(a) ]] )
    return basis

  def get_basis(self, fname):
    acbasis = self.read_ac(fname)
    self.add_caps(acbasis)
    self.acbasis = acbasis

  def HLB20(self, mol):
    auxmol = pyscf.df.make_auxmol(mol, self.acbasis)
    pmol  = mol + auxmol
    eri3c = pmol.intor('int3c2e_sph', shls_slice=(0,mol.nbas,0,mol.nbas,mol.nbas,mol.nbas+auxmol.nbas))
    eri3c = eri3c.reshape(mol.nao_nr(), mol.nao_nr(), -1)
    iao = 0
    for iat in range(auxmol.natm):
      q = auxmol._atom[iat][0]
      for prim in auxmol._basis[q]:
        eri3c[:,:,iao] *= prim[1][1]
        iao+=1
    return numpy.einsum('pqi->pq', eri3c)

  def Heff(self, mol):
    self.mol = mol
    self.Hcore = mol.intor('int1e_nuc_sph') + mol.intor('int1e_kin_sph')
    self.H    = self.Hcore + self.HLB20(mol)
    return self.H

