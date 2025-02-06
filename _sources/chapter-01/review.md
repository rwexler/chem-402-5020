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

### Section 1.1. Course Introduction

- **Thermodynamic Systems**
  - **System vs. Surroundings**
  - **Types**: isolated (no exchange), closed (energy exchange only), open (energy + matter exchange)
- **State Variables and Equilibrium**
  - **State variables** (e.g., $P$, $V$, $T$) vs. **state functions** (e.g., enthalpy $H$, internal energy $U$)
  - **Path functions** (e.g., heat $q$, work $w$) depend on how a process is carried out
  - **Equilibrium** implies no net change in macroscopic properties over time
- **Energy, Heat, and Work**
  - **Kinetic vs. Potential energy**
  - **Heat**: energy transfer due to temperature difference
  - **Work**: energy transfer when a force acts through a distance

### Section 1.2. Kinetic Theory

- **Kinetic Theory Assumptions**
  1. Large number of particles (identical, random motion)
  2. Negligible particle volume compared to container volume
  3. Collisions are elastic (no energy loss)
  4. No long-range interparticle forces (except during collisions)
  5. Newton’s laws (classical mechanics) apply
- **Pressure Origin**
  - Microscopic collisions of particles with container walls
  - Force per unit area results from changes in particle momentum
- **Equipartition Theorem**
  - Average translational kinetic energy per particle $\langle E_{\text{kin}} \rangle = \frac{3}{2} k_{\text{B}} T$
  - Relates temperature to average particle velocity

### Section 1.3. Ideal Gases

- **Classical Gas Laws**
  - **Boyle** ($P \propto 1/V$ at constant $T,\,N$)
  - **Charles** ($V \propto T$ at constant $P,\,N$)
  - **Gay-Lussac** ($P \propto T$ at constant $V,\,N$)
  - **Avogadro** ($V \propto N$ at constant $P,\,T$)
- **Ideal Gas Equation**
  - $PV = N k_{\text{B}} T = n R T$
  - Assumes point-like, non-interacting particles
- **Absolute Temperature Scale**
  - Kelvin scale derived by extrapolating volumes to zero at $-273.15\,^\circ\text{C}$

### Section 1.4. Real Gases

- **Deviations from Ideal Behavior**
  - Become significant at high $P$, high density, or low $T$
  - **Compressibility factor** $Z = \frac{PV}{nRT}$ indicates deviation ($Z=1$ ideal; $Z<1$ net attraction; $Z>1$ net repulsion)
- **van der Waals Equation**
  - Introduces parameters $a$ (accounts for attractions) and $b$ (excluded volume)
  - Shows characteristic “loop” (liquid–vapor coexistence) below critical temperature
- **Critical Point and Corresponding States**
  - Critical conditions $(T_c, P_c, V_{\text{m},c})$ mark the end of the liquid–vapor coexistence line
  - Reduced variables $T_r = T/T_c,\; P_r = P/P_c,\; V_{m,r} = V_m / V_{\text{m},c}$ collapse different substances onto similar curves (principle of corresponding states)

---

## Checklist of Most Important Equations

1. **Ideal Gas Law**

   $$
   PV = N k_{\text{B}} T = n R T.
   $$
   - **Applicability**: low pressures, relatively high temperatures, or low densities (particles effectively non-interacting).

2. **Pressure from Kinetic Theory**

   $$
   P = \frac{N m \langle v^2 \rangle}{3V}.
   $$
   - **Applicability**: idealized gas of point particles with elastic collisions; derived under kinetic-theory assumptions.

3. **Average Kinetic Energy / Equipartition**

   $$
   \langle E_{\text{kin}} \rangle
   \;=\; \frac{3}{2} k_{\text{B}} T
   \quad\Longleftrightarrow\quad
   \frac{1}{2} m \langle v^2 \rangle
   = \frac{3}{2} k_{\text{B}} T.
   $$
   - **Applicability**: classical (high-temperature) regime, ignoring quantum effects and molecular internal modes.

4. **Root-Mean-Square (rms) Speed**

   $$
   v_{\text{rms}}
   = \sqrt{\frac{3 k_{\text{B}} T}{m}}.
   $$
   - **Applicability**: same assumptions as equipartition (kinetic theory of gases).

5. **Compressibility Factor**

   $$
   Z = \frac{PV}{nRT}.
   $$
   - **Interpretation**:
     - $Z = 1$: ideal gas
     - $Z < 1$: net attractive forces
     - $Z > 1$: net repulsive forces
   - **Applicability**: any real gas to quantify deviation from ideality.

6. **van der Waals Equation of State** (in molar form)

   $$
   \left(P + \frac{a_m}{V_m^2}\right)\;\left(V_m - b_m\right)\;=\;R\,T.
   $$
   - **Applicability**: real gases at moderate deviations from ideality; fails under extreme conditions (very high $P$, near liquefaction, etc.).

7. **Critical-Point Relationships (van der Waals)**

   $$
   V_{\text{m},c} = 3\,b_m,\quad
   P_c = \frac{a_m}{27\,b_m^2},\quad
   T_c = \frac{8\,a_m}{27\,b_m\,R}.
   $$
   - **Interpretation**: defines the critical temperature $T_c$, critical pressure $P_c$, and critical molar volume $V_{\text{m},c}$ for a van der Waals fluid.

8. **Corresponding States** (reduced variables)

   $$
   \left(P_r + \frac{3}{V_{m,r}^2}\right)\;\left(V_{m,r} - \frac{1}{3}\right)\;=\;\frac{8}{3}\,T_r
   \quad
   \text{where}
   \quad
   P_r=\frac{P}{P_c},\quad
   T_r=\frac{T}{T_c},\quad
   V_{m,r}=\frac{V_m}{V_{\text{m},c}}.
   $$
   - **Applicability**: near or above critical conditions for fluids that approximate van der Waals behavior; demonstrates “universal” behavior across substances when scaled by critical properties.
