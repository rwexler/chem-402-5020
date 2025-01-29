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

# Ideal Gases

## Overview

Ideal gases are a central model in thermodynamics. By simplifying microscopic interactions and treating molecules as non-interacting point particles, we can accurately predict many macroscopic properties for gases under appropriate conditions. This section reviews the classical gas laws, derives the ideal gas equation of state, and discusses the assumptions underlying the ideal gas model.

---

(section:gas-laws)=

## Gas Laws

```{code-cell} ipython3
:tags: [hide-input]

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import k, N_A
from labellines import labelLines
from myst_nb import glue

# Constants
T = 290  # Temperature in Kelvin (Boyle's Law)
P = 1.0  # Fixed pressure in bar (Charles' Law)
V = 24.53  # Fixed volume in liters (Gay-Lussac's Law)
N = N_A   # Number of particles for Avogadro's Law

temperatures = np.linspace(273.15, 373.15, 101)   # Temperatures in Kelvin
volumes = np.linspace(1, 40, 400)  # Volumes in liters
number_of_particles = np.linspace(0.1, 10, 100) * N_A

# Calculate pressures, volumes, etc., according to each law
pressures_boyle = N * k * T / volumes * 0.01
volumes_charles = (N * k / P * 0.01) * temperatures
pressures_gay_lussac = (N * k / V) * temperatures * 0.01
volumes_avogadro = (k * T / P * 0.01) * number_of_particles

def plot_law(ax, x, y, label, title, xlabel, ylabel):
    line = ax.plot(x, y, "b-", label=label)
    ax.set_title(title, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.grid(True, linestyle="--", linewidth=0.5)
    labelLines(line, zorder=2.5)

fig, axs = plt.subplots(1, 4, figsize=(12, 3))

# Boyle's Law
plot_law(axs[0], volumes, pressures_boyle, "$P = c_\\text{B} / V$", "Boyle (1662)",
         "Volume (L)", "Pressure (bar)")
axs[0].text(0.1, 0.9, 'Constant $T$ & $N$', transform=axs[0].transAxes)

# Charles's Law
plot_law(axs[1], temperatures, volumes_charles, "$V = c_\\text{C} T$", "Charles (1787)",
         "Temperature (K)", "Volume (L)")
axs[1].text(0.1, 0.9, 'Constant $P$ & $N$', transform=axs[1].transAxes)

# Gay-Lussac's Law
plot_law(axs[2], temperatures, pressures_gay_lussac, "$P = c_\\text{GL} T$", "Gay-Lussac (1802)",
         "Temperature (K)", "Pressure (bar)")
axs[2].text(0.1, 0.9, 'Constant $V$ & $N$', transform=axs[2].transAxes)

# Avogadro's Law
plot_law(axs[3], number_of_particles / N_A, volumes_avogadro, "$V = c_\\text{A} N$", "Avogadro (1811)",
         "Number of Particles ($N_\\text{A}$)", "Volume (L)")
axs[3].text(0.1, 0.9, 'Constant $P$ & $T$', transform=axs[3].transAxes)

plt.tight_layout()
glue('gas-laws', fig, display=False)
plt.close(fig)
```

```{glue:figure} gas-laws
:name: fig-gas-laws
:figwidth: 100%
:align: center

Four classical gas laws (blue) shown in chronological order: Boyle's law,[^1] Charles's law, Gay-Lussac's law, and Avogadro's law.

[^1]: [Boyle's original data](https://web.lemoyne.edu/~giunta/classicalcs/boyleform.html)
```

All gases follow these relationships when at least one of the following conditions is met:

- Pressure is relatively low,
- Density is relatively low, and/or
- Temperature is sufficiently high.

These conditions collectively ensure minimal intermolecular interactions, allowing the gas to behave *ideally*.

---

## Deriving an Equation of State

In [Section 1](section-01.md), we defined an **equation of state** as a mathematical relationship among state variables. Above, each gas law relates *two* state variables (among $P, V, T, N$) under conditions where two other variables remain constant. To derive a single equation of state that relates *all* four variables, we can combine these laws and use multivariate calculus (as covered in [Math 233](https://math.wustl.edu/calculus-iii)).

### Total Differential

The total differential $df$ of a function $f$ of $m$ variables $x_1, \ldots, x_m$ is

```{math}
:label: total-differential
df = \sum_{i=1}^m \left( \frac{\partial f}{\partial x_i} \right)_{\{ x_j | j \neq i \}} dx_i,
```

where $\left(\partial f / \partial x_i\right)_{\{ x_j | j \neq i \}}$ is the **partial derivative** of $f$ with respect to $x_i$, holding all other variables constant, and $dx_i$ represents an **infinitesimal change** in $x_i$.

### Total Differential of Volume

From Boyle's, Charles's, and Avogadro's laws, we can treat $V$ as a function of $P, T,$ and $N$:

```{math}
V = V(P, T, N).
```

Applying Equation {eq}`total-differential` to $V$:

```{math}
dV = \left( \frac{\partial V}{\partial P} \right)_{T,N} dP \;+\; \left( \frac{\partial V}{\partial T} \right)_{P,N} dT \;+\; \left( \frac{\partial V}{\partial N} \right)_{P,T} dN.
```

### Partial Derivatives via Gas Laws

Using the gas laws in differential form, one finds:

```{math}
dV \;=\; -\frac{V}{P}\,dP \;+\; \frac{V}{T}\,dT \;+\; \frac{V}{N}\,dN.
```

Rearranging and employing logarithmic differentials:

```{math}
d \ln P \;+\; d \ln V \;=\; d \ln T \;+\; d \ln N.
```

### Integrating the Total Differential

Integrate from an initial state $(P_i, V_i, T_i, N_i)$ to a final state $(P_f, V_f, T_f, N_f)$. By properties of logarithms, we obtain

```{math}
\frac{P_f V_f}{N_f T_f} \;=\; \frac{P_i V_i}{N_i T_i}.
```

Because both initial and final states are arbitrary, the ratio $PV / NT$ must be a constant. We thus arrive at the **ideal gas equation of state**:

```{math}
:label: ideal-gas-equation-of-state
PV \;=\; N k_\text{B} T \;=\; nRT,
```

where $n = N / N_\text{A}$ is the number of moles and $R = k_\text{B} N_\text{A}$ is the [molar gas constant](https://physics.nist.gov/cgi-bin/cuu/Value?r) ($R = 8.314\,\text{J mol}^{-1}\text{K}^{-1}$).

---

## Ideal Gas Assumptions

A gas described by Equation {eq}`ideal-gas-equation-of-state` is called *ideal* because, under low pressures, low densities, or high temperatures, we can adopt simplifying assumptions from kinetic theory:

1. Particles have negligible volume (point particles).
2. Particles experience no intermolecular forces except during elastic collisions.
3. Collisions conserve total energy and momentum.

Many thermodynamic properties derive neatly from these assumptions.

---

## Estimating Particle Distances

Rearranging Equation {eq}`ideal-gas-equation-of-state`, the average distance $\langle d \rangle$ between particles is:

```{math}
:label: average-distance-between-particles
\langle d \rangle = \left(\frac{V}{N}\right)^{1/3} \;=\; \left(\frac{k_\text{B} T}{P}\right)^{1/3}.
```

This simple expression shows that increasing temperature or reducing pressure *increases* the typical separation between gas particles.

````{admonition} Example: Distance Between Air Molecules at 300 K and 1 bar
:class: dropdown

Using Equation {eq}`average-distance-between-particles`:

```{math}
\langle d \rangle 
= \left( \frac{k_\text{B} \cdot 300\,\text{K}}{1.01325 \times 10^5\,\text{Pa}} \right)^{1/3}
\approx 3.45 \times 10^{-9}\,\text{m} \;=\; 34.5\,\text{Å}.
```

In air, Ar atoms (about 0.934% of air molecules) can be modeled with the Lennard-Jones potential:

```{math}
E_\text{pot}(r) 
= 4\epsilon \,\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right],
\quad \epsilon = 0.0103\,\text{eV},\;\; \sigma = 3.4\,\text{Å}.
```

At $r = \langle d \rangle$,

```{math}
E_\text{pot} \approx -3.8\times 10^{-8}\,\text{eV},
```

which is small compared to the thermal kinetic energy,

```{math}
\langle E_\text{kin} \rangle 
\approx \frac{3}{2} k_\text{B} T \;\approx\; 0.039\,\text{eV}.
```

Because $\lvert E_\text{pot} \rvert$ is negligible relative to $\langle E_\text{kin} \rangle$, the Ar atoms behave almost as if they are non-interacting, validating the ideal gas approximation.
````

## Absolute Temperature Scale

```{code-cell} ipython3
:tags: [hide-input]

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.constants import N_A
import statsmodels.api as sm
from myst_nb import glue

# Experimental data for O2 at 1 Pa:
df = pd.read_table("../_static/section-03/isobaric-properties-for-oxygen.tsv", sep="\t").iloc[:, [0, 1, 3]]
df.columns = ["T_C", "P", "V"]
df["v"] = df["V"] * 1e27 / N_A  # Convert volume to nm³/molecule

X = sm.add_constant(df["v"])  # For intercept in linear regression
y = df["T_C"]
model = sm.OLS(y, X).fit()

df["T_C_pred"] = model.predict(X)

fig, ax = plt.subplots(figsize=(4, 4))
ax.plot(df["v"], df["T_C"], "b.")
intercept, slope = model.params
line = ax.plot(df["v"], df["T_C_pred"], "r-", label=f"Slope: {slope:.2e}; Intercept: {intercept:.2f} °C")
labelLines(line, zorder=2.5)

ax.set_xlabel("Volume (nm³/molecule)")
ax.set_ylabel("Temperature (°C)")
ax.grid(True, linestyle="--", linewidth=0.5)

plt.tight_layout()
glue('absolute-temperature', fig, display=False)
plt.close(fig)
```

```{glue:figure} absolute-temperature
:name: fig-absolute-temperature
:figwidth: 100%
:align: center

Linear relationship between the per-molecule volume $v$ and Celsius temperature for O₂ at 1 Pa.[^2]  
The red line is an ordinary least squares fit whose intercept at -273.15 °C signals the theoretical *zero-volume* limit, defining absolute zero (0 K).

[^2]: [NIST WebBook: O₂ Data](https://webbook.nist.gov/cgi/inchi/InChI%3D1S/O2/c1-2)
```

Charles's law states that, at constant (low) pressure, the volume per particle is directly proportional to the gas temperature. Historically, the Celsius scale was used ($0\,^\circ\text{C}$ as the freezing point of water and $100\,^\circ\text{C}$ as the boiling point). However, by measuring gas volume at various Celsius temperatures and extrapolating to where volume vanishes, we identify $-273.15\,^\circ\text{C}$ as the limit—an unphysical negative volume indicates a *theoretical* bound. Shifting the Celsius scale by $273.15$ degrees places absolute zero at $0\,\text{K}$:

```{math}
T(\text{K}) \;=\; T(^\circ\text{C}) \;+\; 273.15.
```

This defines the **Kelvin scale**, an absolute temperature scale that starts at the lowest physically meaningful temperature.
