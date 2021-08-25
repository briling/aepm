import numpy
import pyscf

def readmol(fname, basis, charge=0, spin=0, ignore=False):
  def _makemol(xyz, basis, charge=0, spin=0):
    mol = pyscf.gto.Mole()
    mol.atom = xyz
    mol.charge = charge
    mol.spin = spin
    mol.basis = basis
    mol.build()
    return mol
  def _readxyz(fname):
    with open(fname, "r") as f:
      xyz = f.readlines()
    return "".join(xyz[2:])
  xyz = _readxyz(fname)
  if not ignore:
    mol = _makemol(xyz, basis, charge, spin)
  else:
    try:
      mol = _makemol(xyz, basis)
    except:
      mol = _makemol(xyz, basis, -1)
  return mol

