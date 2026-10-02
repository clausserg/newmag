# OpenMolcas

**$\texttt{NewMag}$** was developed for version 25.06 of OpenMolcas. We do not know whether the program works with other versions (this has not been tested). Therefore, we strongly recommend using **version 25.06** of OpenMolcas.

Regarding OpenMolcas specifically, we strongly recommend placing the files `$Project.rasscf.h5` and `$Project.rassi.h5` in the same directory as the OpenMolcas output file. This ensures more accurate results by utilizing the number of digits available in the `.h5` files.

---

## Minimal Active Space

For a minimal CAS, the OpenMolcas input should be structured as follows:

```console
>>> EXPORT MOLCAS_PRINT = 3 ****REQUIRED****
&SEWARD
RELATIVITY=R02O02 ****REQUIRED****
AMFI ****REQUIRED****
Basis set
...
end of basis
End of input

&SCF
Charge=X
Spin=X
PRORbitals=1 2
End of input

&RASSCF
Spin=X
Symmetry=1 ****REQUIRED****
nActEl=X 0 0
Inactive=XXX
RAS2 = 3, 5 ou 7
CIRoots=X X 1
Iter=200 100
LUMORB
ORBListing=all
ORBAppear=compact
PRWF=0
End of input

&LOCALISATION ****BLOCK REQUIRED FOR THE LOCALISATION AND PURIFICATION OF ACTIVE ORBITALS****
NFrozen=XXX
NORbitals=3, 5 ou 7
CHOLesky
end of input  

&RASSCF
cionly ****REQUIRED****
Spin=X
Symmetry=1 ****REQUIRED**** 
nActEl=X 0 0
Inactive=XXX
RAS2 = 3,5 ou 7
CIRoots=X X 1
Iter=200 100
LUMORB
ORBListing=all ****REQUIRED****
ORBAppear=compact ****REQUIRED****
PRWF=0 ****REQUIRED****
End of input


>>> EXPORT MOLCAS_PRINT = 5 ****REQUIRED****
&RASSI
NROF=1 all
SPINORBIT ****REQUIRED**** 
THRS=0 ****REQUIRED**** 
SOCOupling=0 ****REQUIRED****
end of input
```

All keywords followed by a `****REQUIRED****` comment are **essential** for the **$\texttt{NewMag}$** extraction process.

---
### Key Features of OpenMolcas

OpenMolcas performs two `&RASSCF` blocks:
1. The first computes the CASSCF solutions.
2. An orbital localization is applied to the orbitals of interest, followed by a new CASCI calculation (second `&RASSCF` block).

We encourage users to explore OpenMolcas outputs on GitHub to understand the effect of localization on the wavefunction and CSFs (Configuration State Functions). Initially, the CASSCF yields a highly heterogeneous atomic basis that is difficult to interpret. By leveraging rotational invariance within the active space, orbitals can be localized and purified to express the CSFs in terms of relevant atomic orbitals (p, d, or f).

#### Example Output (Minimal Active Space)
Here is a typical desired result (output of cerocene in the minimal active space, available [here](../output/Molecular/f1/Cerocene/cerocene_CAS_1_7.out)):

**Active orbitals after localization by the `&LOCALISATION` block:**

```console
86    0.0000    0.0000           
       61 CE1    4f1-   ( 0.9989)
87    0.0000    0.0000
       67 CE1    4f1+   ( 0.9989)
88    0.0000    0.0000
       64 CE1    4f0    ( 0.9988)
89    0.0000    0.0000
       55 CE1    4f3-   ( 0.9971)
90    0.0000    0.0000
       73 CE1    4f3+   ( 0.9971)
91    0.0000    0.0000
       58 CE1    4f2-   ( 0.9963)
92    0.0000    0.0000
       70 CE1    4f2+   ( 0.9963)
```

For the first root (for example), the corresponding CSFs are:
```console
printout of CI-coefficients larger than  0.00 for root  1
energy=   -9469.038374
conf/sym  1111111     Coeff  Weight
       1  u000000   0.00000 0.00000
       2  0u00000  -0.00000 0.00000
       3  00u0000   1.00000 1.00000
       4  000u000   0.00000 0.00000
       5  0000u00  -0.00000 0.00000
       6  00000u0   0.00000 0.00000
       7  000000u   0.00000 0.00000
```

Here, the notation `u000000` or `00u0000` represents the active orbitals. In this case, the `ml` order is `-1 +1 0 -3 +1 -2 +2`. Consequently, for this root, the electrons are **only in the $4f_0$ orbital**.

The `ml` order is specified in the OpenMolcas output and is read by **$\texttt{NewMag}$**, though localization may occasionally fail. In such cases, use the `-ml` keyword to manually specify the `ml` order.

> **Note:** For **$\texttt{NewMag}$** extraction, a **CASCI calculation is required**; otherwise, **$\texttt{NewMag}$** cannot identify the wavefunction.

For CASPT2, **$\texttt{NewMag}$** can only read results from a **SA-CASCI/SS-CASPT2/SOCI** output. In other words, the output must contain **only one** `&RASSI` block. An example CASPT2 input is provided below:

```console
...

&RASSCF
cionly ****REQUIRED****
Spin=X
Symmetry=1 ****REQUIRED**** 
nActEl=X 0 0
Inactive=XXX
RAS2 = 3,5 ou 7
CIRoots=X X 1
Iter=200 100
LUMORB
ORBListing=all ****REQUIRED****
ORBAppear=compact ****REQUIRED****
PRWF=0 ****REQUIRED****
End of input

&CASPT2
Multistate=all
MAXIter=100
CONVergence=1.0d-08
PRWF=0
NoMult ****REQUIRED****
End of input

>>> COPY $Project.JobMix JOB001 ****REQUIRED****
>>> EXPORT MOLCAS_PRINT = 5 ****REQUIRED****
&RASSI
EJOB ****REQUIRED****
NROF=1 all
SPINORBIT ****REQUIRED**** 
THRS=0 ****REQUIRED**** 
SOCOupling=0 ****REQUIRED****
end of input
```

> **Important:** For OpenMolcas and the PT2 method, **only the SA-SS-CASPT2 method is implemented** in $\texttt{NewMag}$.

---

## Extended Active Space

For an extended CAS, the principle remains the same as for the minimal CAS, but the `-cas_mo` keyword is **essential**. Use it to select the correct CSFs for the appropriate model space.

### Use Cases for `-cas_mo`
This keyword is particularly useful in scenarios such as:
- Transition $f^{n}d^0 \rightarrow f^{n-1}d^{1}$.
- Including correlation in the active space (with ligand orbitals).

In such cases, CSFs may include orbitals outside the model space. Users must analyze CSF composition to distinguish between model space orbitals and others.

#### Example Output (Extended Active Space)
Here is a typical desired result (output of cerocene in the extended active space, available [here](../output/Molecular/f1/Cerocene/cerocene_CAS_5_9.out)):

**Active orbitals after localization on $4f$ orbitals:**

```console
84    0.0000    1.9998                                                                                              
       58 CE1    4f2-   ( 0.1011) 145 C2     2pz    ( 0.3236) 295 C7     2pz    ( 0.3236) 325 C8     2pz    (-0.3236)
      355 C9     2pz    (-0.3236) 385 C10    2pz    ( 0.3236) 415 C11    2pz    ( 0.3236) 445 C12    2pz    (-0.3236)
      475 C13    2pz    (-0.3236)
85    0.0000    1.9998                                                                                               
       70 CE1    4f2+   ( 0.1011) 115 C1     2pz    (-0.3236) 175 C3     2pz    ( 0.3236) 205 C4     2pz    (-0.3236)
      235 C5     2pz    (-0.3236) 265 C6     2pz    (-0.3236) 505 C14    2pz    ( 0.3236) 535 C15    2pz    ( 0.3236)
      565 C16    2pz    ( 0.3236)
86    0.0000    0.0000                                                                                               
       61 CE1    4f1-   ( 0.9989)
87    0.0000    0.0000                                                                                               
       67 CE1    4f1+   ( 0.9989)
88    0.0000    0.0000                                                                                               
       64 CE1    4f0    ( 0.9989)
89    0.0000    0.0000                                                                                               
       55 CE1    4f3-   ( 0.9970)
90    0.0000    0.0000                                                                                               
       73 CE1    4f3+   ( 0.9970)
91    0.0000    0.0000                                                                                               
       58 CE1    4f2-   ( 0.9961)
92    0.0000    0.0000                                                                                               
       70 CE1    4f2+   ( 0.9961)
```

Here, orbitals 84 and 85 are **$\pi$ orbitals** (based on composition and occupation). Each root has **1890 distinct CSFs** with this orbital distribution: $[\pi \pi f f f f f f f]$. Among these, the CSFs corresponding to the **$f$-orbital model space** must be selected. The seven CSFs to retain for each root are:

```console
1   22u000000
4   220u00000
21  2200u0000
76  22000u000
77  220000u00
78  2200000u0
79  22000000u
```

Here, the goal is to achieve the **same distribution** across relevant orbitals as in the minimal CAS. The `-cas_mo` keyword automates CSF selection but requires arguments specifying the desired organization. In this case:

- The first two orbitals in the CSF always have an occupation of **2**.
- The arguments for `-cas_mo` are `2d 7c`.

  - `2d` sets the leading `2` in the CSFs (e.g., `2d` → `22`).
  - `7c` indicates the **7 $f$-orbitals of interest** (`u000000`, `0u00000`, `00u0000`, etc.).

Combining these arguments automatically selects the **seven relevant CSFs**.

### `-cas_mo` Keyword Categories
The `-cas_mo` keyword supports four categories: `d`, `c`, `o`, and `v` (see the [main README](../README.md) for details):

- **`d`**: Sets a value of `2` in the composition.
- **`c`**: Corresponds to orbitals of interest.
- **`v`**: Assigns the value `0` in the composition (useful for transitions like $f^{1}d^0 \rightarrow f^{0}d^{1}$ or with a virtual orbitals in the active space).
- **`o`**: Assigns values `u`, `d`, `0`, or `2` indiscriminately (useful for transitions like $f^{n}d^0 \rightarrow f^{n-1}d^{1}$).
