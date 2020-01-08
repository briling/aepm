
Programs in this repository work
with basis set files in the Priroda format ([description](https://aip.scitation.org/doi/suppl/10.1063/1.5082231/suppl_file/format.txt))
with some limitations.
In this directory one can find basis sets
used to generate all the examples provided:
* `L1_b2.in` – the L1 set taken from \[[2](../README.md#user-content-ref2)\]
* `L1_b2u.in` – the L1 set taken from \[[3](../README.md#user-content-ref3)\]


If you want to generate starting orbitals with some other basis:
* download it from [Basis Set Exchange](https://www.basissetexchange.org/) in the Dalton format
* run the script, e.g.
  ```
  python3 ./basisconv.py somebasis.dalton > somebasis.in
  ```
* go back and use this basis, e.g.
  ```
  ./q basis/somebasis.in molecule.in
  ```

