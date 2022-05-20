import numpy as np
from pyscf_ext import readmol
from LB2020guess import LB2020guess

molfile = 'HCCCONH2.xyz'
basis   = 'ccpvdz'

mol = readmol(molfile, basis=basis, ignore=True)
Heff = LB2020guess().Heff(mol)

print(Heff)
