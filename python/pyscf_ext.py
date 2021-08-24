import numpy
import pyscf

def makemol(xyz, basis, charge=0, spin=0):
  mol = pyscf.gto.Mole()
  mol.atom = xyz
  mol.charge = charge
  mol.basis = basis
  mol.cart = True
  mol.build()
  return mol

def readmol(fname, basis, charge=0, spin=0, ignore=False):
  def _readxyz(fname):
    with open(fname, "r") as f:
      xyz = f.readlines()
    return "".join(xyz[2:])
  xyz = _readxyz(fname)
  if not ignore:
    mol = makemol(xyz, basis, charge, spin)
  else:
    try:
      mol = makemol(xyz, basis)
    except:
      mol = makemol(xyz, basis, -1)
  return mol

