# Interpretation of $\texttt{NewMag}$ Output

$\texttt{NewMag}$ consists of several distinct sections, which are described in detail in the following subsections.

---

## Structural Analysis

This block calculates:
- The **point group** of the molecule as a whole.
- The **local environment** of the metal center.

It then suggests a **geometric orientation** that allows restarting the CAS calculation using this optimized geometry.

> **Note 1:** $\texttt{NewMag}$ calculates properties based on the XYZ orientation of the molecule specified in the output file (OpenMolcas or ORCA).

> **Note 2:** This block may crash or produce incorrect results (this applies only to this specific block).
---

## Reading Wavefunctions and Energies

This block reads: The **wavefunction** and the **energies** from the output file (OpenMolcas or ORCA)

The selected roots (SR and SO) are displayed according to the calculation levels:
- CASSCF
- CASPT2
- NEVPT2

> **Note :** If your molecular output contains more roots than necessary, the **pseudo-L approximation** is automatically applied to model the ground state.
> Tow alternative approaches are available:
> 
> 1. **Manual root selection**
   Use the keywords `--sr` and/or `--so` to select specific roots of interest. This method applies only to ground state root selection, as our approach cannot model excited states.
   > 
> 1. **Free ion projection method**
   Perform a separate calculation on the free ion using a minimal CAS with the same number of roots as your molecular system (excluding irrelevant roots).
   A projection between the two results will then be carried out to identify the correct ground state roots.

---
## Crystal Field Extraction

This block displays:
- The **effective Hamiltonians** for SR and SO models
- The **composition of the wavefunction** in the Hamiltonian basis

Results are shown according to the calculation levels:
- CASSCF-SR/SO
- CASPT2-SR/SO
- NEVPT2-SR/SO

For SO Hamiltonians, a **symbolic block matrix** is displayed. Below is a general example:

| < $J_{1}$ \| $J_{1}$ ><br>$(1)$ | < $J_{2}$ \| $J_{1}$ ><br>$(j+1)$ | < $J_{3}$ \| $J_{1}$ ><br>$(2j+1)$ | $\ldots$ |
|:----------------------------------:|:----------------------------------:|:----------------------------------:|:------:|
| < $J_{1}$ \| $J_{2}$ ><br>$(2)$ | < $J_{2}$ \| $J_{2}$ ><br>$(j+2)$ | < $J_{3}$ \| $J_{2}$ ><br>$(2j+2)$ | $\ldots$ |
| < $J_{1}$ \| $J_{3}$ ><br>$(3)$ | < $J_{2}$ \| $J_{3}$ ><br>$(j+3)$ | < $J_{3}$ \| $J_{3}$ ><br>$(2j+3)$ | $\ldots$ |
| $\vdots$ | $\vdots$ | $\vdots$ | $\ddots$ |
| < $J_{1}$ \| $J_{j}$ ><br>$(j)$ | < $J_{2}$ \| $J_{j}$ ><br>$(2j)$ | < $J_{3}$ \| $J_{j}$ ><br>$(3j)$ | $\ldots$ |

The blocks are then displayed in groups of columns, in the order indicated by the brackets in the matrix:
1. (1)
2. (2)
3. (3)
...

> **Note 1:** Hamiltonian values are expressed in $\text{cm}^{-1}$.
> 
> **Note 2:** Wavefunction composition values are expressed as percentages.
---
## Extraction of Crystal Field Parameters

This block:
1. Calculates all CFP parameters composing each Hamiltonian using the ITO method
2. Reads/extracts and converts parameters if the output file contains:
   - `AILFT` elements
   - `single_aniso` elements

A Markdown table is generated listing the $B_k^q$ values according to each calculation level. Below are all calculation levels with their corresponding extraction programs:

|         | casscf-ailft     | casscf-so         | casscf-sr         | nevpt2-ailft     | nevpt2-sr         | nevpt2-so         | caspt2-sr         | caspt2-so         | x_aniso_J                                | x_aniso_L                                |
| :-----: | :--------------- | :---------------- | :---------------- | :--------------- | :---------------- | :---------------- | :---------------- | :---------------- | :--------------------------------------- | :--------------------------------------- |
| program | $\texttt{AILFT}$ | $\texttt{NewMag}$ | $\texttt{NewMag}$ | $\texttt{AILFT}$ | $\texttt{NewMag}$ | $\texttt{NewMag}$ | $\texttt{NewMag}$ | $\texttt{NewMag}$ | $\texttt{Single-aniso}$ | $\texttt{Single-aniso}$ |

> **Note 1:** Values are expressed in cm⁻¹ and, following Stevens convention, are unnormalized.
> 
> **Note 2:** Values $S_k$, $S^q$ and $S$ are calculated according to Stevens convention and are normalized.
> 
> **Note 3:**
> - `x_aniso`: Refers to the xth iteration of the $\texttt{Single}\_\texttt{aniso}$  program
> - For $\texttt{AILFT}$: Only the first iteration shown in the output is considered
---
## Model Spectra

This block reconstructs the energy spectrum (SR and SO) using the $B_k^q$ values from the previous block.

For the SO level, the ITO method is used to extract the The SOC constant $\lambda$ and $\zeta$ (spherical approximation)
