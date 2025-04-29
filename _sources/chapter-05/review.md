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

1. **Origins & Roles of Free Energies**
   - **Helmholtz Free Energy**
     Defined as

     ```{math}
     A \;=\; U - T\,S,
     ```

     with natural variables $\left(T,V\right)$. At fixed temperature and volume, a process is spontaneous if $\Delta A<0$.
   - **Gibbs Free Energy**
     Defined as

     ```{math}
     G \;=\; H - T\,S \;=\; U + P\,V - T\,S,
     ```

     with natural variables $\left(T,P\right)$. At constant $T,P$, spontaneity requires $\Delta G<0$, and $\Delta G=0$ marks equilibrium.
   - **Statistical‐Mechanical Link**
     The partition function $Q$ underpins both; e.g.,

     ```{math}
     A = -k_B T \ln Q,
     ```

     connecting macroscopic free energy to microscopic states.

2. **The Third Law of Thermodynamics**
   - **Planck’s Statement**:
     As $T\to0$, the entropy of any **pure, defect-free crystalline** substance approaches zero.
   - **Boltzmann’s View**:
     $S = k_B\ln\Omega$, where $\Omega$ is the count of microstates. A perfect crystal has $\Omega=1$ at 0 K.
   - **Residual Entropy**:
     Real crystals often contain equivalent configurations (e.g., orientational disorder in solid CO), giving a nonzero “residual” $S(0)$. Example: $S_\text{res}=k_B\ln2$ per molecule.

3. **T- and P-Dependence of $\Delta G$**
   - **Temperature Effects**:

     ```{math}
     \Delta G(T_f)-\Delta G(T_i)
     = \int_{T_i}^{T_f}\!\Delta C_p\,dT \;-\; T\,\int_{T_i}^{T_f}\!\frac{\Delta C_p}{T}\,dT.
     ```

     Splits enthalpic and entropic contributions; used to assess ammonia formation from 298 K to 773 K, showing $\Delta G$ becomes less negative (and eventually positive) at high $T$.
   - **Pressure Effects**:
     At fixed $T$,

     ```{math}
     dG = V\,dP,
     ```

     so for an ideal‐gas reaction with $\Delta n$ change in moles,

     ```{math}
     \Delta G(P) = \Delta G^\circ + RT\,\Delta n\,\ln\frac{P}{P^\circ}.
     ```

     Example: raising $P$ to ∼319 bar at 773 K makes ammonia synthesis spontaneous again.

---

## 2. Checklist of Most Important Equations

| Equation | Variables | Narrative & Use Cases | Caveats/Approximations |
|:-|:-|:-|:-|
| $dU = T\,dS - P\,dV + \sum_i \mu_i\,dN_i$ | $U$: internal energy<br>$S$: entropy<br>$P$: pressure<br>$V$: volume<br>$\mu_i$, $N_i$: chemical potentials & particle numbers | Fundamental First Law in natural variables. Basis for deriving all thermodynamic potentials. | Omits non‐PV work; assumes simple homogeneous system. |
| $H = U + P\,V,\quad dH = T\,dS + V\,dP$ | $H$: enthalpy | Convenient for constant‐$P$ processes (e.g., reactions in open vessels). | Requires only PV work. |
| $A = U - T\,S,\quad dA = -S\,dT - P\,dV$ | $A$: Helmholtz free energy | Governs spontaneity under fixed $T,V$. Maximum non‐PV work equals $-\Delta A$. | Valid for closed systems at constant $T,V$. |
| $G = H - T\,S,\quad dG = -S\,dT + V\,dP$ | $G$: Gibbs free energy | Determines whether processes at constant $T,P$ are spontaneous ($\Delta G<0$). | Holds for constant‐$T,P$; neglects non‐PV work. |
| $A = -k_B T \ln Q$ | $k_B$: Boltzmann constant<br>$Q$: partition function | Links macroscopic $A$ to microscopic energy levels; compute equilibrium constants and thermodynamic functions.  | Exact only if full $Q$ known; often approximated. |
| $S = \frac{U}{T} + k_B\ln Q$ | — | Statistical‐mechanical entropy in the canonical ensemble.  | Same as above. |
| $S = k_B\ln \Omega$ | $\Omega$: number of microstates | Boltzmann’s principle for isolated systems; underlies the Third Law.  | Conceptual; $\Omega$ intractable in practice. |
| $\Delta S = \displaystyle\int_{0}^{T} \frac{C_p(T')}{T'}\,dT'$ | $C_p$: heat capacity at constant pressure | From Third Law: computes $S$ from 0 K (where $S=0$ for a perfect crystal).  | Assumes no phase transitions; $C_p\to0$ as $T\to0$. |
| $\Delta H = \displaystyle\int_{T_i}^{T_f} C_p(T)\,dT$ | — | Temperature corrections to standard‐state enthalpies.  | Needs $C_p(T)$; often treated as constant. |
| $dG = -S\,dT + V\,dP$ | — | Basis for Clapeyron relation and pressure/temperature dependence of equilibrium. | Requires known $S$ and $V$. |
| $G(P) = G^\circ + RT\ln\frac{P}{P^\circ}$ | $R$: gas constant | Adjusts Gibbs energy for nonstandard pressures in ideal gas approximation; used to find equilibrium pressure.  | Ideal‐gas assumption. |
| $\displaystyle\Delta G(P,T) = \Delta G^\circ(T) + RT\,\Delta n\,\ln\frac{P}{P^\circ}$ | $\Delta n$: net change in moles of gas | Solves for $P$ or $T$ at which a gas‐phase reaction becomes spontaneous or reaches equilibrium.  | Assumes all species ideal gases; single pressure for all components. |
