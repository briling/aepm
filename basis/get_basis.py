import pyscf

basis = 'minao'
qs = [1, 6, 8, 9, 15, 17, 27, 45, 77]

print(f'''$basis
type=general
default={basis}
set={basis}''')
for q in qs:
    bas = pyscf.gto.basis.load(basis, pyscf.data.elements.ELEMENTS[q])
    lenb = sum([len(shell[1])-1 for shell in bas])
    print("      %3d %3d"%(q, lenb))
    for shell in bas:
        l = shell[0]
        for i in range(1, len(shell[1])):
            print(l, len(shell)-1)
            for prim in shell[1:]:
                print(prim[0], prim[i])
print('$end')
