# ORCA

$\texttt{NewMag}$ was developed for **version 6.1.0** of ORCA. While we cannot guarantee compatibility with other versions (as this has not been tested), we strongly recommend using **version 6.1.0** of ORCA for optimal results.

---

## Minimal Active Space

For a minimal CAS calculation in ORCA, the input should be structured as follows:

```console
! nori 
! autoaux 

%rel
  method DKH 
  picturechange 2
end

%casscf
  nel     X
  norb    X
  mult    X
  nroots  X
  maxiter 120 
  printwf det   ####REQUIRED####
  actorbs forbs ####REQUIRED#### (or: actorbs dorbs)
  ci            
   TPrintWF   0 ####REQUIRED####
  end           
  rel           
    dosoc true  ####REQUIRED####
    tprint 0.0  ####REQUIRED####
  end
end

*xyz X X
GEOM_XYZ
*
```

All keywords followed by a `####REQUIRED####` comment are **essential** for the $\texttt{NewMag}$ extraction process.

---
### Key Features of ORCA

The `actorbs` keyword automatically localizes and purifies the orbitals of interest at the end of the CASSCF calculation. The generated orbital set is standardized with the following $m_l$ orbital order:

- For `forbs`: `0 +1 -1 +2 -2 +3 -3`
- For `dorbs`: `0 +1 -1 +2 -2`
#### Example Output (Minimal Active Space)
Here is a typical desired result (output of europium complex in minimal active space, available [here](../output/Molecular/f6/)):

**Active orbitals after the CASSCF calculation:**

```console
                    186       187       188       189       190       191   
                  -0.35127  -0.34602  -0.31888   0.00619   0.00611   0.00857
                   2.00000   2.00000   2.00000   0.85714   0.85714   0.85714
                  --------  --------  --------  --------  --------  --------
 0 Eu f0              0.0       0.0       0.0      99.8       0.0       0.0
 0 Eu f+1             0.0       0.0       0.0       0.0      99.8       0.0
 0 Eu f-1             0.0       0.0       0.0       0.0       0.0      99.5
 2 O  px              0.3       0.0      13.3       0.0       0.0       0.0
 3 O  px              0.1      13.2       0.4       0.0       0.0       0.0
 4 O  px              5.7       0.2       3.9       0.0       0.0       0.0
20 N  px              6.5       0.0       1.8       0.0       0.0       0.0
22 C  px              0.2       0.2      13.3       0.0       0.0       0.0
23 C  px             11.7       0.5      11.7       0.0       0.0       0.0
24 C  px             24.7       0.2       0.0       0.0       0.0       0.0
25 C  px              0.9       0.2      20.3       0.0       0.0       0.0
26 C  px             10.9       0.5       7.8       0.0       0.0       0.0
27 C  px             19.9       0.0       8.6       0.0       0.0       0.0
33 C  px              0.0       7.0       0.2       0.0       0.0       0.0
34 C  px              0.2       9.6       0.1       0.0       0.0       0.0
35 C  px              0.2      18.4       0.2       0.0       0.0       0.0
37 C  px              0.2      12.5       0.1       0.0       0.0       0.0
38 C  px              0.2      15.0       0.1       0.0       0.0       0.0

                    192       193       194       195       196       197   
                   0.00693   0.00598   0.00559   0.00601   0.03076   0.05908
                   0.85714   0.85714   0.85714   0.85714   0.00000   0.00000
                  --------  --------  --------  --------  --------  --------
 0 Eu f+2            99.7       0.0       0.0       0.0       0.0       0.0
 0 Eu f-2             0.0      99.9       0.0       0.0       0.0       0.0
 0 Eu f+3             0.0       0.0      99.8       0.0       0.0       0.0
 0 Eu f-3             0.0       0.0       0.0      99.7       0.0       0.0
20 N  px              0.0       0.0       0.0       0.0       0.1      12.4
21 N  px              0.0       0.0       0.0       0.0      14.0       0.4
22 C  px              0.0       0.0       0.0       0.0       0.2       9.9
24 C  px              0.0       0.0       0.0       0.0       0.2      15.7
26 C  px              0.0       0.0       0.0       0.0       0.1       5.4
28 C  px              0.0       0.0       0.0       0.0       0.5      20.9
32 C  px              0.0       0.0       0.0       0.0      22.5       0.3
33 C  px              0.0       0.0       0.0       0.0       5.4       0.2
34 C  px              0.0       0.0       0.0       0.0       6.8       0.2  
36 C  px              0.0       0.0       0.0       0.0      15.9       0.5
38 C  px              0.0       0.0       0.0       0.0      10.1       0.2
```

For the first root (for example), the corresponding CSFs are:
```console
ROOT   0:  E=  -14754.8277991795 Eh
   [uuuuuu0]     -0.102826281
   [uuuuu0u]     -0.017272320
   [uuuu0uu]     -0.178089435
   [uuu0uuu]     -0.189160472
   [uu0uuuu]     -0.950904872
   [u0uuuuu]     -0.088049553
   [0uuuuuu]     -0.098275591
```

Here, the notation `uuuuuu0` or `uuuuu0u` represents the active orbitals. In this case, the $m_l$ order is `0 +1 -1 +2 -2 +3 -3`. Consequently, for this root the electrons are mainly in the orbitals: $4f_{0}$, $4f_{+1}$, $4f_{+2}$, $4f_{-2}$, $4f_{+3}$ and $4f_{-3}$

> **Important:** For ORCA and the PT2 method, **only the NEVPT2 method is implemented** in $\texttt{NewMag}$.

---
## Extended Active Space

For an extended CAS with ORCA, the approach remains similar to OpenMolcas (see [documentation](./OpenMolcas_input.md)). However, the method for obtaining this type of active space using the `actorbs` keyword differs.

We recommend performing three consecutive calculations:

1. **Initial CASSCF Calculation**
   Perform a CASSCF calculation in your desired final active space
   Output: `wfn_mix.gbw`

2. **Orbital Localization**
   Using the previous wavefunction file (`wfn_mix.gbw`), perform a CASCI calculation:
   - In the minimal active space
   - With the `actorbs` keyword
   - With `maxiter 1`
   Output: `wfn_loc_CAS_minimal.gbw`

3. **Final CASCI Calculation**
   Using the localized wavefunction file (`wfn_loc_CAS_minimal.gbw`), perform another CASCI calculation:
   - In your desired final active space
   - With `maxiter 1`
   Output: Localized and purified orbitals in extended active space

After completing these steps, you should obtain an extended active space containing the localized and purified orbitals of interest.
