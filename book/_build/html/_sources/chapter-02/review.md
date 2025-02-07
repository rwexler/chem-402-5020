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

## Checklist of Key Concepts

### Section 2.1: Introduction to Statistical Mechanics

1. **Microscopic–Macroscopic Connection**  
   - Macroscopic (thermodynamic) properties can be understood as **statistical averages** of microscopic properties.
   - “Expected value” (ensemble average) is the central idea: $\langle X \rangle = \sum X_i\,p_i$.

2. **Arithmetic Average vs. Expected Value**  
   - Arithmetic average: $\bar{X} = \frac{1}{M}\sum X_i$.  
   - Expected value: $\langle X \rangle = \sum X_i\,p_i$, where $p_i$ is the probability of microstate $i$.

3. **Microstates and Ensembles**  
   - **Microcanonical ensemble**: $(N, V, E)$ — all accessible microstates have the same $E$.  
   - **Canonical ensemble**: $(N, V, T)$ — system can exchange energy with a heat bath.  
   - **Grand canonical ensemble**: $(\mu, V, T)$ — system exchanges both energy and particles with the reservoir.

4. **Fundamental Postulate (Microcanonical)**  
   - For an isolated system, **each accessible microstate is equally probable**.

### Section 2.2: Canonical Ensemble

1. **Closed System**  
   - Exchanges energy (heat) with surroundings; cannot exchange matter.

2. **Boltzmann Factor and Partition Function**  
   - Probability of microstate $i$:  

     $$
       p_i = \frac{e^{-\beta E_i}}{Q},
       \quad
       \beta = \frac{1}{k_{\text{B}} T}.
     $$  
   - **Partition function** $Q$: normalizes probabilities; sum of Boltzmann factors over all states:

     $$
       Q = \sum_{i} e^{-\beta E_i}.
     $$

3. **Two-State System**  
   - Simplest nontrivial example: energies $E_1$ and $E_2$.  
   - Partition function: $Q = e^{-\beta E_1} + e^{-\beta E_2}$.  
   - Probabilities:
   
     $$
       p_1 = \frac{1}{1 + e^{-\beta\,\Delta E}},
       \quad
       p_2 = 1 - p_1,
       \quad
       \Delta E = E_2 - E_1.
     $$

4. **Interpretation of $Q$**  
   - “Effective number” of thermally accessible microstates; as $T \to 0$, few states are accessible; as $T \to \infty$, many states are accessible.

### Section 2.3: Ensemble Averages

1. **Internal Energy**  
   - $U \equiv \langle E \rangle$.  
   - In canonical ensemble:

     $$
       U = \frac{1}{Q}\sum_i E_i\, e^{-\beta E_i} = -\left(\frac{\partial \ln Q}{\partial \beta}\right)_{N,V}.
     $$

2. **Heat Capacity ($C_V$)**  
   - $C_V = \bigl(\frac{\partial U}{\partial T}\bigr)_{N,V}$.  
   - Also relates to **energy fluctuations**:
   
     $$
       \sigma_E^2 \;=\; \langle(E - \langle E\rangle)^2\rangle \;=\; k_{\text{B}}\,T^2\,C_V.
     $$

3. **Pressure** (brief introduction)  
   - Later formal definition involves $\bigl(\frac{\partial U}{\partial V}\bigr)_{S}$.  
   - In canonical ensemble, $\displaystyle P \;=\; k_{\text{B}}\,T\,\bigl(\frac{\partial\ln Q}{\partial V}\bigr)_{N,T}$.

---

## 2. Checklist of Most Important Equations

1. **Expected Value of a General Variable**  

   $$
     \langle X \rangle \;=\;\sum_i X_i\,p_i
     \quad\text{or}\quad
     \int X(\omega)\,p(\omega)\,d\omega.
   $$  
   - **Applicability**: general definition in statistical mechanics/probability theory.

2. **Microcanonical Ensemble (Fundamental Postulate)**  

   $$
     p_i \;=\;\frac{1}{M}
     \quad\text{(for all accessible microstates)}.
   $$  
   - **Applicability**: isolated system with fixed $(N, V, E)$.

3. **Canonical Ensemble Probability**  

   $$
     p_i \;=\;\frac{e^{-\beta\,E_i}}{Q},
     \quad
     \beta = \frac{1}{k_{\text{B}}\,T},
   $$  
   - **Applicability**: closed system in thermal equilibrium with a heat bath at temperature $T$.

4. **Canonical Partition Function**  

   $$
     Q \;=\;\sum_i e^{-\beta\,E_i}.
   $$  
   - **Applicability**: all systems that can be described by $(N, V, T)$; sum or integral runs over all microstates.

5. **Internal Energy in the Canonical Ensemble**  

   $$
     U \;=\;\langle E\rangle
     \;=\;\frac{1}{Q}\,\sum_i E_i\,e^{-\beta E_i}
     \;=\;-\left(\frac{\partial \ln Q}{\partial \beta}\right)_{N,V}.
   $$  
   - **Applicability**: same as above (canonical ensemble).

6. **Heat Capacity at Constant Volume**  

   $$
     C_V
     \;=\;
     \left(\frac{\partial U}{\partial T}\right)_{N,V}.
   $$  
   - **Applicability**: measures how internal energy changes with $T$ when $N, V$ are constant.  

7. **Energy Fluctuations**  

   $$
     \sigma_E^2
     \;=\;\langle E^2\rangle - \langle E\rangle^2
     \;=\;
     k_{\text{B}}\,T^2\,C_V.
   $$  
   - **Applicability**: canonical ensemble; relates variance in energy to the heat capacity.

8. **Pressure (Canonical)**  

   $$
     P
     \;=\;
     k_{\text{B}}\,T \,\left(\frac{\partial \ln Q}{\partial V}\right)_{N,T}.
   $$  
   - **Applicability**: canonical ensemble; a precursor to more rigorous derivation using thermodynamic potentials.
