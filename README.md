
Codes for the paper
**"Atomic effective potentials for starting molecular electronic structure calculations"** \[[1]\]

## Contents

* [Build](#build)
* [Usage](#usage)
   * [1. Generate starting orbitals](1-generate-starting-orbitals)
   * [2. Optimize caps' exponents](2-optimize-caps-exponents)
* [Files](#files)
* [References](#references)


## Build

```
make q qap
```
### Requirements:
* `GNU/Linux` or `Cygwin`
* `gcc >= 4.7`


## Usage

### 1. Generate starting orbitals

To generate starting orbitals for a molecule, run
```
./q <basis> <molecule>.{in,out} [options]
```
| command-line option      |  description                                                                                                                                                       | default value                                                               |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------- |
| `vectors:%s`             | file to save the vectors in                                                                                                                                  | `<molecule>.vec`                                                                  |
| `print:%d`               | `1`: print density matrix; `2`: print density matrix in the lower-triangular form; `3`: print vectors in column-major order; `4`: print vectors in `MOLDEN` format | `0`                                                                         |
| `aaar:%d`                | if > 0, use the two-component scalar-relativistic approximation \[[3]\], the parameters have to be specified in the basis set file                                 | `1` if the parameters are specified                                         |
| `finite_nuclei:%d`       | if > 0, use the finite Gaussian nucleus model \[[4]\]                                                                                                              | `1` if the scalar-relativistic approximation is enabled; <br> `0` otherwise |

#### Examples

to be found in [mol/](mol/):
```
./q basis/L1_b2.in  mol/H5C5FeC5H5.{in,out} print:2
./q basis/L1_b2.in  mol/Cr2N2O2C8O8.{in,out} print:2
./q basis/L1_b2.in  mol/H6C4FeC3O3.{in,out} print:2
./q basis/L1_b2u.in mol/H8C8ThC8H8.{in,out} print:2
./q basis/L1_b2.in  mol/TEMPO.{in,out} print:2
./q basis/L1_b2.in  mol/HCCCONH2.{in,out} print:4
```

#### Vectors file format

```
The binary file contains:
  1   int32_t --- number of basis functions M
  M   doubles --- orbital energies (alpha)
  M*M doubles --- orbital vectors  (alpha)
  M   doubles --- orbital energies (beta)
  M*M doubles --- orbital vectors  (beta)
Thus its size is sizeof(int32_t)+2*M*(M+1)*sizeof(double) bytes.

Spherical basis functions are used and their order is as follows:
  P: P0, P+1, P-1
  D: D0, D+1, D-1, D+2, D-2
  F: F0, F+1, F-1, F+2, F-2, F+3, F-3
etc.
```

### 2. Optimize caps' exponents
To optimize the
<img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;a_{n+1}" title="a_{n+1}" height="14" />
values on a set of molecules, run
```
./qap <basis> <inputfile> [outputfile] [options]
```
| command-line option      |  description                                               |                                                                                                                                                                                                                                                                                   default value(s)         |
| ------------------------ | ---------------------------------------------------------- |                                                                                                                                                                                                                                                                                   ------------------------ |
| `f:%s`                   | measure to use for optimization: `S`, `E`, `S0`, `E0`      |                                                                                                                                                                                                                                                                                   `S`                      |
| `check:%d`               | if > 0, instead of optimization, compute the S and E measures on the molecules: <br> `1` – use HF-based parameters and caps with the default or specified in the input file values of <img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;a_{n+1}" title="a_{n+1}" height="14" />; `2` – use HFS-based parameters |    `0`                      |
| `testgrad:%d`            | if > 0, compare analytical and numerical derivatives of all the measures with respect to the <img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;a_{n+1}" title="a_{n+1}" height="14" /> values, and quit |                                                                                                            `0`                      |
| `save:%d`                | if > 0, save guess vectors after optimization              |                                                                                                                                                                                                                                                                                   `0`                      |
| `np:%d`                  | number of threads                                          |                                                                                                                                                                                                                                                                                   `1`                      |
| `o:%d,%lf`               | optimization parameters (maximum number of iterations and convergence criterion) |                                                                                                                                                                                                                                                             `128,1e-5`               |
| `o1:%d,%lf,%lf,%lf,%lf`  | various 1D optimization parameters (maximum number of iterations, convergence criterion, default value of second derivative estimation, maximum step length, minimum value of second derivative estimation) |                                                                                                                                  `32,1e-6,3.0,0.125,1e-3` |

#### `qap` input file format

See [qap_ex/](qap_ex) for examples.

First, the initial values of
<img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;a_{n+1}" title="a_{n+1}" height="14" />
can be specified:
```
value1  element number(s)
value2  element number(s)
...
valueK  element number(s)
```
For example,
```
    0.333   1-2
    0.0625  3
    0.0625  4
    0.333   5-10
fix 0.333   13
    0.125   14-18
```
sets
<img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;a_{n+1}=0.333" title="a_{n+1}=0.333" height="14" />
for Hydrogen, Helium, and Boron through Neon;
<img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;0.0625" title="0.0625" height="14" />
for Lithium and Berillium;
<img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;0.125" title="0.125" height="14" />
for Aluminum through Argon; and the default value of
<img src="http://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;0.125" title="0.125" height="14" />
for all other elements.
During the optimization, the parameters within the groups (in this case, `1-2`, `5-10`, `11-12`, `14-18`, and `19-102`)
are considered to be equal.
One can freeze some of the parameters by starting lines with the `fix` keyword.

Next, molecular geometry files are listed (comments starting with `#` are allowed).
For example,
```
H2.in
H2O.in # water
/directory/with/molecules/NH3.in
# try later:
# Se8.in
S8.in
```
File names should correspond the mask `%.in`
and for each one a vector file called `%.vec` should exist.

#### Examples

to be found in [mol_opt/](mol_opt/) and [qap_ex/](qap_ex):
```
./qap basis/L1_b2.in  qap_ex/light.{in,opt.out}  np:4  # optimize cap exponents on 4 processors

./qap basis/L1_b2.in  qap_ex/light.{in,f.out} check:1  # compute S and E measures

./qap basis/L1_b2u.in qap_ex/SeTePo.{in,out} np:4      # optimize cap exponents for Se,Te,Po
                                                       # with fixed exponents for H,O,S
```


## Files

* `src/`            – source directory
* `obj/`            – build directory
* `mol/*.in`        – molecular geometry files
* `mol/*.vec`       – starting vectors computed with `q`
* `mol/*.out`       – `q` output files
* `mol_opt/*/*.in`  – molecular geometry files
* `mol_opt/*/*.vec` – corresponding SCF vectors computed with Priroda-19 (PBE/L1)
* `qap_ex/`         – `qap` input and output files
* `basis/`          – basis sets taken from \[[2]\] and \[[3]\]


## References

<a name="ref1">\[1\]</a>
D. N. Laikov and K. R. Briling, [Theor. Chem. Acc.][LB2020] **139**, 17 (2020).

<a name="ref2">\[2\]</a>
D. N. Laikov, [Theor. Chem. Acc.][L2019b] **138**, 40 (2019).

<a name="ref3">\[3\]</a>
D. N. Laikov, [J. Chem. Phys.][L2019a] **150**, 061103 (2019).

<a name="ref4">\[4\]</a>
L. Visscher and K. G. Dyall, [At. Data Nucl. Data Tables][VD1997] **67**, 207 (1997)

[1]: #user-content-ref1
[2]: #user-content-ref2
[3]: #user-content-ref3
[4]: #user-content-ref4
[LB2020]:https://doi.org/10.1007/s00214-019-2521-3
[L2019b]:https://doi.org/10.1007/s00214-019-2432-3
[L2019a]:https://doi.org/10.1063/1.5082231
[VD1997]:https://doi.org/10.1006/adnd.1997.0751

