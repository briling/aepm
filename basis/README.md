
Programs in this repository work
with basis set files in the Priroda format ([description](https://aip.scitation.org/doi/suppl/10.1063/1.5082231/suppl_file/format.txt))
with some limitations.
In this directory one can find basis sets
used to generate all the examples provided:
* `L1_b2.in` – the L1 set taken from \[[2](../README.md#user-content-ref2)\]
* `L1_b2u.in` – the L1 set taken from \[[3](../README.md#user-content-ref3)\]

If you want to generate starting orbitals for some other basis set:
* download it from [Basis Set Exchange](https://www.basissetexchange.org/) in the Dalton format and
* run the conversion [script](basisconv.py).

For example, to use STO-3G one can run in the main directory:
  ```
  python3 basis/basisconv.py <((curl "https://www.basissetexchange.org/api/basis/sto-3g/format/dalton/?version=1&optimize_general=true")) > basis/sto3g.in
  ./q basis/sto3g.in your-molecule.in
  ```

Another way is to use the `get_basis.py` script (if you have PySCF).
