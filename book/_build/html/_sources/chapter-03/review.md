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

### Section 3.1: Conservation of Energy

1. **Overview**
   - Energy conservation underlies all of thermodynamics.
   - Early demonstrations (e.g., Bernoulli, Euler) merged Newton’s mechanics with the idea that total mechanical energy remains constant in an isolated system.

2. **Mechanical Equivalent of Heat**
   - Joule’s experiments showed that mechanical work can be converted into heat in a fixed, quantitative ratio.
   - For water: $\text{MEH} \approx 4.18\,\mathrm{J}\,\mathrm{g}^{-1}\,^\circ\mathrm{C}^{-1}$, establishing heat as another form of energy transfer.

3. **First Law of Thermodynamics**
   - Statement: $\Delta U = q + w$, where $q$ is heat **absorbed by** the system, $w$ is work **done on** the system.
   - **Path vs. State Functions**: $U$ is a state function (depends only on the current state), while $q$ and $w$ are path functions (depend on how the process is carried out).
   - **Sign Convention**: in chemistry, $w>0$ if work is done **on** the system.

4. **Types of Work**
   - Generalized form: $\delta w = \text{(generalized force)} \times \text{(generalized displacement)}$.
   - Common examples in thermodynamics:
     - **$P\,dV$ work**: compression/expansion of a gas ($\delta w = -P\,dV$).
     - **Surface tension**: $\delta w = -\gamma\,dA$.
     - **Hooke’s law**: $\delta w = k\,(x - x_0)\,dx$.
   - Example calculations:
     - Expansion against constant pressure.
     - Inflating a bubble with surface tension.
     - Stretching a spring-like fiber obeying Hooke’s law.

---

### Section 3.2: Applications of the First Law

1. **Thermodynamic Processes**
   - **Quasi-static**: carried out slowly enough for the system to remain in (near) equilibrium at each step.
   - **Reversible**: idealized, quasi-static + no dissipative losses; can be undone with no net change in system or surroundings.
   - **Irreversible**: any real process that breaks one or more reversibility conditions.

2. **How to Apply the First Law**
   - **Stepwise Procedure**:
     1. Choose two independent variables (e.g., $V, T$ or $P, T$).
     2. Rewrite $\delta q$ and $\delta w$ in terms of differentials of those variables.
     3. Apply relevant constraints (e.g., isobaric, isochoric, isothermal, adiabatic).
     4. Use the system’s equation of state (e.g., $PV=nRT$ for an ideal gas).
     5. Integrate or sum changes to find total $q$, $w$, and $\Delta U$.

3. **Using $V$ and $T$ as Independent Variables**
   - **First Law**: $\delta q = dU + P\,dV$.
   - For an ideal gas ($\partial U/\partial V)_T = 0$):

     ```{math}
       \delta q 
       = \left(\frac{\partial U}{\partial T}\right)_V dT + P\,dV.
     ```

   - Special processes:
     - **Isochoric** ($dV=0$): $q = \int C_V\,dT$.
     - **Isothermal** ($dT=0$): $q = \int P\,dV = nRT \ln(V_2/V_1)$.
     - **Adiabatic** ($q=0$): $\Delta U = w$. For a monoatomic ideal gas, $T\,V^{\tfrac{2}{3}}=\text{const}$.

4. **Microscopic Interpretation**
   - **Internal Energy**: $U = \sum_i p_i E_i$.
   - **Heat ($q$)** changes probabilities $p_i$.
   - **Work ($w$)** shifts the energy levels $E_i$ themselves (e.g., compressing the container changes spacing of quantum states).
   - This perspective unifies macroscopic thermodynamics with microscopic statistical mechanics.

---

### Section 3.3: Enthalpy

1. **Definition of Enthalpy**
   - $H = U + PV$.
   - For a process at constant pressure, $\delta q_P = dH$.

2. **$\Delta H$ and Heat at Constant Pressure**
   - $\Delta H = q_P$. Calorimetry at $P=\text{const}$ measures enthalpy changes directly.
   - **Heat Capacity at Constant Pressure**: $C_P = \bigl(\tfrac{\partial H}{\partial T}\bigr)_P$.

3. **Standard Enthalpies**
   - **Standard State**: $P^\circ=1\text{ bar}$; tabulated data often at $298.15\,\mathrm{K}$.
   - **Standard Enthalpy of Formation**, $\Delta H_f^\circ$: enthalpy change to form 1 mole of a compound from its elements in their standard states.
     - By convention, $\Delta H_f^\circ=0$ for any element in its standard state.
   - **Standard Enthalpy of Reaction**, $\Delta H_{\mathrm{rxn}}^\circ$:

     ```{math}
       \Delta H_{\mathrm{rxn}}^\circ
       \;=\;\sum_{\text{products}} \nu_p H_p^\circ 
       \;-\;\sum_{\text{reactants}} \nu_r H_r^\circ.
     ```

   - **Hess’s Law**: enthalpy changes are path independent; $\Delta H$ of an overall reaction is the algebraic sum of enthalpy changes for its steps.

---

## 2. Checklist of Most Important Equations

1. **First Law of Thermodynamics**

   ```{math}
   \Delta U \;=\; q \;+\; w,
   ```

   - $q$: heat absorbed by the system.
   - $w$: work done on the system.
   - **Differential form**: $dU = \delta q + \delta w.$

2. **Pressure–Volume Work (Constant External $P$)**

   ```{math}
   w \;=\; -\,P_{\text{ext}}\;\bigl(V_f - V_i\bigr).
   ```

   - Sign convention: expansion ($\Delta V>0$) $\Rightarrow$ $w<0$ (system does work on surroundings).

3. **Heat for Isochoric Process**

   ```{math}
   q_{V} \;=\; \int C_V \, dT
   \;\;\Longrightarrow\;\;
   q_{V} = C_V\,\Delta T
   \;\;\text{if }C_V\approx\text{const.}
   ```

4. **Heat for Isothermal Process (Ideal Gas)**

   ```{math}
   q_{\text{isoT}} 
   \;=\; \int_{V_i}^{V_f} P\,dV
   \;=\; n\,R\,T\,\ln{\bigl(\tfrac{V_f}{V_i}\bigr)}.
   ```

5. **Adiabatic Condition (No Heat Exchange)**

   ```{math}
   q_{\text{adi}} = 0
   \;\;\Longrightarrow\;\;
   \Delta U = w.
   ```

   - For a monoatomic ideal gas: $T\,V^{\,2/3} = \text{constant}.$

6. **Enthalpy Definition**

   ```{math}
   H \;=\; U + P\,V.
   ```

   - At constant $P$: $\Delta H = q_P.$

7. **Standard Enthalpy of Reaction**

   ```{math}
   \Delta H_{\mathrm{rxn}}^\circ
   \;=\;
   \sum_{p} \nu_p\,H_p^\circ
   \;-\;
   \sum_{r} \nu_r\,H_r^\circ.
   ```

   - **Hess’s Law** and **enthalpies of formation** provide a straightforward way to compute $\Delta H_{\mathrm{rxn}}^\circ$.
