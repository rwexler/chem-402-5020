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

### Section 2.1: Introduction to Statistical Mechanics

1. **Microscopic–Macroscopic Connection**
   - Macroscopic (thermodynamic) properties can be understood as **statistical averages** of microscopic properties.
   - “Expected value” (ensemble average) is the central idea:

   ```{math}
   \langle X \rangle \;=\; \sum_i X_i\,p_i.
   ```

2. **Arithmetic Average vs. Expected Value**
   - Arithmetic average:

   ```{math}
   \bar{X} \;=\; \frac{1}{M}\sum_{i=1}^M X_i.
   ```

   - Expected value:

   ```{math}
   \langle X \rangle \;=\; \sum_i X_i\,p_i,
   \quad p_i = \text{probability of microstate } i.
   ```

3. **Microstates and Ensembles**
   - **Microcanonical ensemble** $(N, V, E)$: system is isolated, fixed total energy $E$.
   - **Canonical ensemble** $(N, V, T)$: system in thermal contact with a reservoir at temperature $T$.
   - **Grand canonical ensemble** $(\mu, V, T)$: system can exchange both energy and particles with a reservoir.

4. **Fundamental Postulate (Microcanonical)**
   - For an isolated system, **each accessible microstate is equally probable**.

---

### Section 2.2: Canonical Ensemble

1. **Closed System**
   - Exchanges energy (heat) with surroundings; no exchange of matter.

2. **Boltzmann Factor and Partition Function**
   - Probability of microstate $i$:

   ```{math}
   p_i \;=\; \frac{e^{-\beta E_i}}{Q},
   \quad \beta \;=\;\frac{1}{k_{\text{B}}\,T}.
   ```

   - **Partition function** $Q$:

   ```{math}
   Q \;=\; \sum_{i} e^{-\beta E_i},
   ```

   which normalizes probabilities.

3. **Two-State System**
   - A simple example with energies $E_1$ and $E_2$.
   - Partition function:

   ```{math}
   Q \;=\; e^{-\beta E_1} + e^{-\beta E_2}.
   ```

   - Probabilities (if $\Delta E = E_2 - E_1$):

   ```{math}
   p_1 = \frac{1}{1 + e^{-\beta\,\Delta E}},
   \quad
   p_2 = 1 - p_1.
   ```

4. **Interpretation of $Q$**
   - $Q$ is like an “effective count” of accessible microstates.
   - At low $T$, only the lowest energy states matter; at high $T$, many states are accessible.

---

### Section 2.3: Ensemble Averages

1. **Internal Energy**
   - $U \equiv \langle E \rangle$.
   - In the canonical ensemble:

   ```{math}
   U
   \;=\;
   \frac{1}{Q} \sum_i E_i\, e^{-\beta E_i}
   \;=\;
   -\,\bigl(\tfrac{\partial \ln Q}{\partial \beta}\bigr)_{N,V}.
   ```

2. **Heat Capacity $(C_V)$**
   - Measures how internal energy changes with temperature:

   ```{math}
   C_V = \bigl(\tfrac{\partial U}{\partial T}\bigr)_{N,V}.
   ```

   - Also related to **energy fluctuations**:

   ```{math}
   \sigma_E^2
   \;=\;
   \langle (E - \langle E\rangle)^2\rangle
   \;=\;
   k_{\text{B}}\,T^2\,C_V.
   ```

3. **Pressure** (brief introduction)
   - In the canonical ensemble:

   ```{math}
   P
   \;=\;
   k_{\text{B}}\,T\,
   \bigl(\tfrac{\partial\ln Q}{\partial V}\bigr)_{N,T}.
   ```

---

### Section 2.4: Molecular Partition Functions

1. **Many Identical and Independent Particles**
   - For $N$ **distinguishable** particles, each in a one‐particle partition function $q$, the total partition function is

   ```{math}
   Q \;=\; q^N
   \quad (\text{if particles are distinct}).
   ```

   - For **indistinguishable** particles (e.g., identical atoms/molecules obeying single occupancy),

   ```{math}
   Q \;=\; \frac{q^N}{N!}.
   ```

2. **Molecular Partition Function**
   - Within Born–Oppenheimer, a molecule’s energy decomposes into translational, rotational, vibrational, and electronic parts:

   ```{math}
   \varepsilon = \varepsilon_{\text{trans}} + \varepsilon_{\text{rot}} + \varepsilon_{\text{vib}} + \varepsilon_{\text{elec}}.
   ```

   - Thus,

   ```{math}
   q
   \;=\;
   q_{\text{trans}}\,
   q_{\text{rot}}\,
   q_{\text{vib}}\,
   q_{\text{elec}}.
   ```

---

### Section 2.5: Particle in a Box

1. **Quantum Levels**
   - For a 1D box $(0 \le x \le L)$,

   ```{math}
   E_n = \frac{h^2}{8mL^2}\,n^2,\;\; n=1,2,3,\dots
   ```

   - In 3D (a rectangular box of sides $L_x, L_y, L_z$):

   ```{math}
   E_{n_x,n_y,n_z}
   \;=\;
   \frac{h^2}{8m}
   \left(\frac{n_x^2}{L_x^2} + \frac{n_y^2}{L_y^2} + \frac{n_z^2}{L_z^2}\right).
   ```

2. **Partition Function (3D Box)**
   - For a cubic box of volume $V = L^3$:

   ```{math}
   q
   \;=\;
   \sum_{n_x=1}^\infty \sum_{n_y=1}^\infty \sum_{n_z=1}^\infty
   \exp\!\Bigl[
    -\beta \frac{h^2}{8mL^2}\,(n_x^2 + n_y^2 + n_z^2)
    \Bigr].
    ```

   - At high $T$ or large $L$, we approximate the sum by an integral and obtain the **classical** partition function:

   ```{math}
   q_{\text{trans}}
   \;=\;
   \frac{V}{\Lambda^3},
   \quad
   \Lambda
   \;=\;
   \sqrt{\frac{h^2}{2\pi m k_{\text{B}} T}}\,\,
   \text{(thermal de Broglie wavelength).}
   ```

3. **Many-Particle System (Ideal Gas)**
   - For $N$ indistinguishable, non-interacting particles:

   ```{math}
   Q
   \;=\;
   \frac{q^N}{N!}
   \;=\;
   \frac{1}{N!}
   \biggl(\frac{V}{\Lambda^3}\biggr)^{\!N}.
   ```

   - Leads directly to the ideal gas law and to:

   ```{math}
   U \;=\; \frac{3}{2}\,N\,k_{\text{B}}T,
   \quad
   C_V \;=\; \frac{3}{2}\,N\,k_{\text{B}},
   \quad
   P = \frac{N\,k_{\text{B}}\,T}{V}.
   ```

---

### Section 2.6: Harmonic Oscillator

1. **Quantum Harmonic Oscillator**
   - Energy levels for a 1D harmonic oscillator of frequency $\omega$:

   ```{math}
   E_n
   \;=\;
   \hbar\omega\,\Bigl(n + \tfrac{1}{2}\Bigr),
   \quad
   n=0,1,2,\dots
   ```

2. **Partition Function**
   - Summing over these energy levels:

   ```{math}
   q
   \;=\;
   \sum_{n=0}^{\infty} e^{-\beta E_n}
   \;=\;
   \frac{e^{-\frac{1}{2}\,\beta\,\hbar\omega}}{1 - e^{-\beta\,\hbar\omega}}.
   ```

3. **Ensemble Averages**
   - **Internal Energy**:

   ```{math}
   U
   \;=\;
   \frac{\hbar\omega}{2}
   \;+\;
   \frac{\hbar\omega}{\,e^{\beta\,\hbar\omega}-1\,}.
   ```

   - **Heat Capacity**:

   ```{math}
   C_V
   \;=\;
   k_{\text{B}}
   \Bigl(\frac{\hbar\omega}{k_{\text{B}} T}\Bigr)^{2}\,
   \frac{\,e^{\hbar\omega/(k_{\text{B}}T)}\,}
        {\bigl(e^{\hbar\omega/(k_{\text{B}}T)} - 1\bigr)^{2}}.
   ```

   - At high temperature ($k_{\text{B}}T \gg \hbar\omega$), each harmonic oscillator recovers the **classical limit** $U \to k_{\text{B}}T$ and $C_V \to k_{\text{B}}$.

<!-- 4. **Einstein Model** (for solids)
   - Considers each atom as a 3D harmonic oscillator with the same $\omega$.
   - Provides a simple route to the Dulong–Petit law ($C_V \approx 3 N k_{\text{B}}$) at high $T$. -->

---

### Section 2.7: Linear Rigid Rotor

1. **Energy Levels**
   - For a rigid, linear rotor of moment of inertia $I$:

   ```{math}
   E_J
   \;=\;
   \frac{\hbar^2}{2I}\,J\,(J+1),
   \quad
   J=0,1,2,\dots
   ```

   - Each level $E_J$ has degeneracy $g_J = 2J + 1$.

2. **Rotational Partition Function**
   - Exact form:

   ```{math}
   q_{\text{rot}}
   \;=\;
   \sum_{J=0}^\infty
   (2J + 1)\, e^{-\beta \,\frac{\hbar^2}{2I}\,J\,(J+1)}.
   ```

   - Define the **rotational temperature** $\Theta_{\text{rot}} = \frac{\hbar^2}{2k_{\text{B}}I}$.
   - **High-temperature limit** ($k_{\text{B}}T \gg \frac{\hbar^2}{2I}$) gives

   ```{math}
   q_{\text{rot}}
   \;\approx\;
   \frac{T}{\Theta_{\text{rot}}}
   \quad (\text{for a heteronuclear diatomic}).
   ```

   - For a **homonuclear** diatomic, include a symmetry factor $\sigma=2$, so

   ```{math}
   q_{\text{rot}}
   \;\approx\;
   \frac{T}{\sigma\,\Theta_{\text{rot}}}.
   ```

3. **Ensemble Averages** (High-$T$ Approximation)
   - **Internal Energy** ($U_{\text{rot}}$):

   ```{math}
   U_{\text{rot}} \;\approx\; k_{\text{B}}\,T.
   ```

     (One linear rotor has 2 rotational degrees of freedom $\rightarrow \frac{2}{2} k_{\text{B}}T$.)
   - **Heat Capacity**:

   ```{math}
   C_V^{(\text{rot})} \;\approx\; k_{\text{B}}.
   ```

---

### Section 2.8: Molecular Statistical Mechanics

1. **Combining All Degrees of Freedom**
   - For a general molecule, the total one‐molecule partition function is

   ```{math}
   q
   \;=\;
   q_{\text{trans}}\,
   q_{\text{rot}}\,
   q_{\text{vib}}\,
   q_{\text{elec}}.
   ```

   - Extends to polyatomic molecules with more complex rotational constants $\Theta_{\text{rot,A}}, \Theta_{\text{rot,B}}, \Theta_{\text{rot,C}}$ and multiple vibrational frequencies $\Theta_{\text{vib},j}$.

2. **Symmetry Considerations**
   - Symmetry factor $\sigma$ must be included for molecules with indistinguishable orientations (e.g., homonuclear diatomics, symmetrical polyatomics).

3. **Summary Table**
   - Often, we tabulate $q_{\text{trans}}, q_{\text{rot}}, q_{\text{vib}}, q_{\text{elec}}$ for different molecule types (linear, nonlinear, spherical top, symmetric top, etc.), applying high-temperature (classical) or more exact quantum results as needed.

---

## 2. Checklist of Most Important Equations

Below is a **unified list** of the major equations from Sections 2.1–2.8.

A. **Expected Value of a General Variable**

```{math}
\langle X \rangle
\;=\;
\sum_i X_i\,p_i
\quad\text{or}\quad
\int X(\omega)\,p(\omega)\,d\omega.
```

- **Applicability**: general definition in statistical mechanics/probability theory.

---

B. **Microcanonical Ensemble (Fundamental Postulate)**

```{math}
p_i
\;=\;
\frac{1}{M}
\quad(\text{for all accessible microstates}).
```

- **Applicability**: isolated system, $(N, V, E)$.

---

C. **Canonical Ensemble Probability**

```{math}
p_i
\;=\;
\frac{\,e^{-\beta\,E_i}\,}{Q},
\quad
\beta = \frac{1}{k_{\text{B}}\,T}.
```

- **Applicability**: closed system in thermal contact at $T$.

---

D. **Canonical Partition Function**

```{math}
Q
\;=\;
\sum_i e^{-\beta\,E_i}.
```

- **Applicability**: $(N, V, T)$ ensemble; sum/integral over all microstates.

---

E. **Internal Energy (Canonical)**

```{math}
U
\;=\;
\langle E\rangle
\;=\;
\frac{1}{Q}\,\sum_i E_i\,e^{-\beta E_i}
\;=\;
-\bigl(\tfrac{\partial \ln Q}{\partial \beta}\bigr)_{N,V}.
```

---

F. **Heat Capacity at Constant Volume**

```{math}
C_V
\;=\;
\bigl(\tfrac{\partial U}{\partial T}\bigr)_{N,V}.
```

Also,

```{math}
\sigma_E^2
\;=\;
k_{\text{B}}\,T^2\,C_V.
```

---

G. **Pressure (Canonical)**

```{math}
P
\;=\;
k_{\text{B}}\,T \,\bigl(\tfrac{\partial \ln Q}{\partial V}\bigr)_{N,T}.
```

---

H. **Molecular Systems: Indistinguishability Factor**

```{math}
Q
\;=\;
\frac{q^N}{N!},
```

- **Applicability**: identical, indistinguishable particles (e.g., ideal gases).

---

G. **Translational Partition Function** (Particle in a 3D box at high $T$)

```{math}
q_{\text{trans}}
\;=\;
\frac{V}{\Lambda^3},
\quad
\Lambda
=
\sqrt{\frac{h^2}{2\pi m k_{\text{B}}T}}.
```

---

I. **Harmonic Oscillator Partition Function (1D)**

```{math}
q_{\text{HO}}
\;=\;
\frac{\,e^{-\tfrac{1}{2}\beta\hbar\omega}\,}{1 - e^{-\beta\,\hbar\omega}}.
```

- **Internal Energy**:

```{math}
U_{\text{HO}}
\;=\;
\frac{\hbar\omega}{2}
+
\frac{\hbar\omega}{\,e^{\beta\,\hbar\omega}-1\,}.
```

---

J. **Rotational Partition Function (Linear Rotor, High $T$)**

```{math}
q_{\text{rot}}
\;\approx\;
\frac{T}{\sigma\,\Theta_{\text{rot}}},
\quad
\Theta_{\text{rot}} = \frac{\hbar^2}{2k_{\text{B}}\,I}.
```

- **Internal Energy**:

```{math}
U_{\text{rot}} \;\approx\; k_{\text{B}}\,T.
```

- **Heat Capacity**:

```{math}
C_V^{(\text{rot})} \;\approx\; k_{\text{B}}.
```

---

K. **Polyatomic Molecules**

- **General form** (neglecting interactions):

```{math}
q
\;=\;
q_{\text{trans}}\;
q_{\text{rot}}\;
q_{\text{vib}}\;
q_{\text{elec}}.
```

- For linear vs. nonlinear rotors or multiple vibrational modes, each factor is included appropriately (with possible symmetry factors).
