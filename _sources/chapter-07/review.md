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

1. **Gibbs Free Energy & Reaction Extent**
   - Define the **extent of reaction** $x$ so that for  
     $\ce{\nu_A A(g) + \nu_B B <=> \nu_Y Y + \nu_Z Z}$, $n_{\ce{A}} = n_{\ce{A},0}-\nu_{\ce{A}} x$, etc.
   - The **driving force** is

     ```{math}
     \left. \frac{\partial G}{\partial x} \right|_{T,P}
     = \sum_i \nu_i\,\mu_i = \Delta G_{\text{rxn}}.
     ```

     - If $\Delta G_{\text{rxn}}<0$, reaction proceeds forward.
     - If $\Delta G_{\text{rxn}}=0$, equilibrium.
     - If $\Delta G_{\text{rxn}}>0$, reaction moves backward.

2. **Reaction Quotient & Equilibrium Constant**
   - **Reaction quotient** $Q$ compares partial pressures to standard state:

     ```{math}
     Q = \prod_i \left( \frac{P_i}{P^{\circ}} \right)^{\nu_i},
     ```

     so $\Delta G_{\text{rxn}} = \Delta G_{\text{rxn}}^\circ + RT\ln Q$.
   - At **equilibrium**, $\Delta G_{\text{rxn}}=0\implies Q=K$, the **equilibrium constant**:

     ```{math}
     K = \exp\left(-\frac{\Delta G_{\text{rxn}}^\circ}{RT}\right).
     ```

   - **Interpretation**: $Q<K$ drives products; $Q>K$ drives reactants.

3. **Gas‐Phase Example & ICE Analysis**
   - Dimerization $\ce{2NO2(g) <=> N2O4(g)}$: use an ICE table with extent $x$ to express mole numbers and mole fractions; define

     ```{math}
     K_p = \frac{P_{\ce{N2O4}}/P^\circ}{(P_{\ce{NO2}}/P^\circ)^2}.
     ```

   - Solve for $x_\text{eq}$ at various temperatures (e.g., 250 K, 298 K, 350 K) and pressures.
   - Observe that increasing $T$ (endothermic dimerization) shifts equilibrium toward reactants.

4. **Temperature Dependence: van't Hoff Equation**
   - From $\Delta G_{\text{rxn}}^\circ = -RT\ln K$ and Gibbs–Helmholtz,

     ```{math}
     \left. \frac{d\ln K}{dT}\right|_P \;=\;\frac{\\Delta H_{\text{rxn}}^\circ}{R\,T^2}.
     ```

   - Enables prediction of how $K$ changes with $T$, and rationalizes the data trend in the $\ce{NO2}$–$\ce{N2O4}$ system.

5. **Microscopic Basis via Partition Functions**
   - Helmholtz energy $A=-RT\ln Q$ for a mixture of ideal gases with $\;Q(T,V,\{N_i\})=\prod_i Q_i(T,V,N_i)$.
   - Chemical potentials $\mu_i=-RT\,(\partial\ln Q/\partial N_i)$, leading to the same equilibrium condition $\sum\nu_i\,\mu_i=0$.
   - From these, one can compute $K$ purely from molecular partition functions, e.g., for $\ce{H2(g) + I2(g) <=> 2HI(g)}$.

---

## 2. Checklist of Most Important Equations

| **Equation** | **Variables** | **Meaning & Use Cases** | **Caveats/Approximations** |
|:-|:-|:-|:-|
| $\Delta G_{\text{rxn}} = \sum_i \nu_i\,\mu_i \;=\;\left(\frac{\partial G}{\partial x}\right)_{T,P}$ | $\nu_i$: stoichiometric coefficients (products positive)<br>$\mu_i$: chemical potential of species $i$ | Reaction driving force. Use this to assess spontaneity and locate equilibrium $\left( \Delta G_{\text{rxn}} = 0 \right)$. | Assumes constant $T,P$ and only $PV$ work; applies generally to closed, simple systems. |
| $\mu_i = \mu_i^\circ + RT\ln\left( P_i / P^{\circ} \right)$ | — | Defines chemical potential in terms of measurable quantities. Controls mass exchange; at equilibrium between phases or compartments, $\mu_i$ equal in each. | Requires proper choice of standard state and correction for nonideal behavior. |
| $\Delta G_{\text{rxn}} = \Delta G_{\text{rxn}}^\circ + RT\ln Q$ | $\Delta G_{\text{rxn}}^\circ$: standard‐state free energy change<br>$Q=\prod \left( P_i / P^{\circ} \right)^{\nu_i}$: reaction quotient | Relates instantaneous driving force to composition. Use to predict direction of reaction under non‐standard conditions. | — |
| $K = \exp\left(-\frac{\Delta G_{\text{rxn}}^\circ}{RT}\right)$ | $K$: equilibrium constant (can be $K_P$, $K_c$, etc.) | Defines equilibrium composition. Use standard thermodynamic data ($\Delta G_{\text{rxn}}^\circ$) to compute $K$. | Assumes temperature‐independent $\Delta G_{\text{rxn}}^\circ$ unless corrected; different $K$ for different standard states (pressure vs concentration). |
| $Q = \frac{\left( y_{\ce{Y}} P / P^{\circ} \right)^{\nu_{\ce{Y}}} \, \left( y_{\ce{Z}} P / P^{\circ} \right)^{\nu_{\ce{Z}}}}{\left( y_{\ce{A}} P / P^{\circ} \right)^{\nu_{\ce{A}}} \, \left( y_{\ce{B}} P / P^{\circ} \right)^{\nu_{\ce{B}}}}$ | $y_i$: mole fraction of gas $i$<br>$P$: total pressure | Specific form of $Q$ for ideal‐gas mixtures. Use in ICE analyses to solve for equilibrium extent $x$. | Ideal‐gas assumption; neglects nonidealities. |
| $\frac{d\ln K}{dT} = \frac{\Delta H_{\text{rxn}}^\circ}{R\,T^2}$ | $\Delta H_{\text{rxn}}^\circ$: standard enthalpy change | **van't Hoff**: predicts how $K$ shifts with $T$. | Assumes $\Delta H_{\text{rxn}}^\circ$ constant over the temperature range; in reality, $C_p$ corrections may be needed. |
| $A = -RT\ln Q,\quad \mu_i = -RT\left(\frac{\partial\ln Q}{\partial N_i}\right)_{T,V}$ | $A$: Helmholtz free energy of the mixture<br>$Q$: total partition function | Connects molecular partition functions to macroscopic free energies and chemical potentials. Enables **microscopic** calculation of $K$. | Exact only if full partition functions (translational, rotational, vibrational, electronic) are known. |
