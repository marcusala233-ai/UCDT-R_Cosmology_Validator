

# UCDT-R* Cosmology: Unified Dark Sector Dynamics

[![DOI](https://zenodo.org/badge/1134625974.svg)](https://doi.org/10.5281/zenodo.18290985)

![Status](https://img.shields.io/badge/Status-Validated-success) ![Physics](https://img.shields.io/badge/Physics-BSM_Cosmology-blueviolet) ![License](https://img.shields.io/badge/License-MIT-green)

> **Abstract:** This repository contains the numerical validation of the **UCDT-R* (Unified Geometric Unification - Refactored)** theory. We demonstrate that a single pNGB scalar field, embedded in a stabilized 5D manifold with signature $(3,2)$, successfully reproduces the phenomenology of both Dark Energy and Cold Dark Matter (CDM). Through cycle-averaged equation of state analysis, we confirm a transition from a repulsive inflationary phase ($\langle w \rangle \approx -1$) to a clustering matter phase ($\langle w \rangle \approx 0$).

![UCDT-R Visual Compendium](ucdt_visual_compendium.png)

## 1. Theoretical Framework

The UCDT-R* model addresses the Dark Sector problem by introducing a pseudo-Nambu-Goldstone Boson (pNGB) $\phi$ coupled to a 5D geometry. To ensure unitarity and avoid Ostrogradsky instabilities typically associated with multi-time signatures, we impose the **Stabilized Non-Minimal Compactification (SNMC)** condition on the energy-momentum tensor:

$$T_{44}^{(5)} \ge 0 \quad (\forall t)$$

The scalar field dynamics are governed by a periodic potential protected by shift symmetry:

$$V(\phi) = \Lambda_{UV}^4 \left[ 1 - \cos\left(\frac{\phi}{f_a}\right) \right]$$

## 2. Methodology: Cycle Averaging Analysis

Direct observation of scalar field dynamics reveals high-frequency oscillations in the Equation of State ($w(t)$). To determine the effective macroscopic behavior of the fluid, we implement a **Cycle Averaging** algorithm.

Instead of instantaneous measurement, we calculate the effective equation of state $\langle w \rangle_{eff}$ over an oscillation period $\tau \sim m_{eff}^{-1}$:

$$\langle w \rangle_{eff} = \frac{\int_{t}^{t+\tau} P(t') dt'}{\int_{t}^{t+\tau} \rho(t') dt'}$$

This method filters out kinetic noise and reveals the true thermodynamic nature of the UCDT-R* fluid.

## 3. Numerical Results

The simulation, performed via an RK4 solver with the parameters $\Lambda_{UV}=1.5$ and $f_a=0.5$, yields the following cosmological evolution:

### Phase I: Dark Energy Dominated ($t < 1.0$)
The field resides in the slow-roll regime at the potential hill.
* **Observed $\langle w \rangle$:** `-0.7665`
* **Physical Interpretation:** Repulsive pressure driving cosmic acceleration (Inflation/DE).

### Phase II: Phase Transition ($1.0 < t < 2.5$)
The field rolls down towards the potential minimum.
* **Observed $\langle w \rangle$:** `+0.25` (Transient Kinetic Regime)
* **Physical Interpretation:** "Reheating" or kinetic-dominated transition phase.

### Phase III: Dark Matter Dominated ($t > 3.0$)
The field oscillates at the bottom of the potential.
* **Observed $\langle w \rangle$:** Converges to `0.00` (e.g., `0.08`, `-0.02`, `-0.11`, `+0.05`).
* **Physical Interpretation:** The pressure averages to zero. The fluid behaves effectively as **Cold Dark Matter (CDM)**, allowing for gravitational collapse and structure formation.

#### Data Log Excerpt (Cycle Averaged):
```text
Tempo(fim)   <w> (EOS)         Regime Efetivo (Macro)
-------------------------------------------------------------
0.78         -0.7665           [ENERGIA ESCURA] (Repulsivo)
...
3.18          0.0858           ** MATERIA ESCURA ** (Aglomeravel)
3.98         -0.0236           ** MATERIA ESCURA ** (Aglomeravel)
9.58         -0.0003           ** MATERIA ESCURA ** (Aglomeravel)
19.98         0.0845           ** MATERIA ESCURA ** (Aglomeravel)

4. Conclusion
The numerical evidence confirms that UCDT-R* provides a consistent unification mechanism. The vanishing effective sound speed ($c_{s,eff}^2 \approx 0$) in the late universe resolves the structure formation issues typically affecting scalar field cosmologies.

🤝 Citação
Se você utilizar este software ou a teoria UCDT-R em sua pesquisa, por favor, cite utilizando o DOI arquivado no Zenodo (botão acima) ou o arquivo CITATION.cff presente neste repositório.

📝 Licença
Distribuído sob a licença MIT. Veja o arquivo LICENSE para mais informações.
Autor: Marcus Ala Pedreira Roriz
Pesquisador em Física Teórica e Cosmologia Computacional
Principal Investigator: Marcus Validation Engine: C++ Custom Solver (RK4 + Cycle Analyzer)
