---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Review

## 1. Checklist of Key Concepts

1. **Phase‐Equilibrium Condition**
   - **Gibbs Equality**: At equilibrium between any two phases $\alpha$ and $\beta$, their molar Gibbs free energies are equal:

     ```{math}
     G_\alpha = G_\beta.
     ```

     This condition holds at a phase‐transition point (e.g. melting point, boiling point) and underpins all phase‐diagram constructions.

2. **Latent Heat & Entropy Jump**
   - **Enthalpy Change**: $\Delta H_\text{trs}$ is nonzero
     <!-- at a first‐order transition -->

     ```{math}
     \Delta H_\text{trs} = H_\beta - H_\alpha,
     ```

     and known as the *latent heat*.
   - **Entropy Change**: $\Delta S_\text{trs}$ likewise jumps:

     ```{math}
     \Delta S_\text{trs} = S_\beta - S_\alpha = \frac{\Delta H_\text{trs}}{T_\text{trs}}.
     ```

3. **Chemical Potential in Open Systems**
   - **Definition**: For component $i$ in an open, multi‐component system,

     ```{math}
     \mu_i \;=\;\Bigl(\frac{\partial G}{\partial N_i}\Bigr)_{T,P,N_{j\neq i}},
     ```

     appears in the fundamental relation

     ```{math}
     dG = -S\,dT + V\,dP + \sum_i \mu_i\,dN_i.
     ```

   - **Role**: $\mu_i$ is the “driving force” for matter exchange across a boundary; at phase equilibrium, the chemical potentials of each component are equal in all phases.

4. **Clausius–Clapeyron Equation**
   - **Differential Form**: Relates the slope of a phase boundary on a $P$–$T$ diagram to latent heat and volume change:

     ```{math}
     \frac{dP}{dT}
     \;=\;
     \frac{\Delta S_\text{trs}}{\Delta V_\text{trs}}
     \;=\;
     \frac{\Delta H_\text{trs}}{T\,\Delta V_\text{trs}}.
     ```

   - **Sketching Phase Diagrams**: By knowing $\Delta H_\text{trs}$ and $\Delta V_\text{trs}$, one can approximate the coexistence curve between two phases (e.g. solid–liquid, liquid–gas).
   - **Integrated Form (Ideal Gas Approximation)**:

     ```{math}
     \ln\frac{P_2}{P_1}
     \approx
     -\frac{\Delta H_\text{vap}}{R}\,\Bigl(\frac{1}{T_2}-\frac{1}{T_1}\Bigr),
     ```

     useful for extracting $\Delta H_\text{vap}$ from vapor‐pressure data.

---

## 2. Checklist of Most Important Equations

<!-- Assumes both phases are in stable equilibrium;  -->
<!-- ; neglects P‐dependence of $H$ -->
<!-- Assumes transition occurs reversibly at a single $T$. -->
<!-- ; can break down near critical point where $\Delta V\to0$ -->
<!-- and phases are well‐mixed -->
<!-- interactions beyond ideal solution or  -->
<!-- Use in heterogeneous phase equilibria (e.g., gas–solid adsorption). -->
<!-- ; use fugacity for real gases at high pressure -->

| **Equation** | **Variables** | **Meaning & Use Cases** | **Caveats/Approximations** |
|:-|:-|:-|:-|
| $G_\alpha = G_\beta$ | $G$: molar Gibbs free energy | **Fundamental** phase‐equilibrium criterion. Use to locate transition points in $P$–$T$ space. | Applies only at the equilibrium curve. |
| $\Delta H_\text{trs} = H_\beta - H_\alpha$ | $H$: enthalpy | **Latent heat** associated with phase change (fusion, vaporization). Determines heat required/​released at transition. | Measured at constant $P$. |
| $\Delta S_\text{trs} = S_\beta - S_\alpha = \dfrac{\Delta H_\text{trs}}{T_\text{trs}}$ | $S$: entropy; $T_\text{trs}$: transition temperature | Quantifies **entropy jump** at a first‐order transition. Essential for Clapeyron equation. | — |
| $\dfrac{dP}{dT} = \dfrac{\Delta S_\text{trs}}{\Delta V_\text{trs}} = \dfrac{\Delta H_\text{trs}}{T\,\Delta V_\text{trs}}$ | $\Delta V_\text{trs} = V_\beta - V_\alpha$: molar volume change | **Clausius–Clapeyron**: slope of coexistence curve. Use to predict how $P$ changes with $T$ for phase boundaries (e.g., vapor pressure curve). | Requires knowledge of $\Delta H_\text{trs}$ and $\Delta V_\text{trs}$. |
| $\ln\frac{P_2}{P_1} \approx -\frac{\Delta H_\text{vap}}{R}\Bigl(\frac{1}{T_2}-\frac{1}{T_1}\Bigr)$ | $R$: gas constant | **Integrated** Clapeyron under ideal‐gas assumption for vapor: use vapor‐pressure measurements at two $T$ to find $\Delta H_\text{vap}$ or predict $P$. | Assumes vapor behaves ideally and $\Delta H_\text{vap}$ is constant over $T_1$–$T_2$. |
| $\mu_i = \Bigl(\frac{\partial G}{\partial N_i}\Bigr)_{T,P,N_{j\neq i}}$ | $\mu_i$: chemical potential of component $i$<br>$N_i$: number of moles | Defines the intensive variable controlling **mass exchange**. Equal in all phases at multi‐phase equilibrium (e.g., for water in coexisting liquid & vapor). | Valid for systems where composition can vary; neglects real gas corrections. |
| $dG = -S\,dT + V\,dP + \sum_i \mu_i\,dN_i$ | — | **Fundamental differential** for open, multi‐component systems. Basis for deriving all equilibrium conditions (thermal, mechanical, chemical). | Assumes only $PV$ work and matter exchange; neglects electric, magnetic, surface‐tension work, etc. |
| $G(P,T) = G^\circ(T) + RT\ln\frac{P}{P^\circ}$ | $P^\circ$: standard pressure | Adjusts the Gibbs free energy of an ideal gas from standard pressure to $P$. | Ideal‐gas assumption. |
