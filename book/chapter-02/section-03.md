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

# 2.3. Ensemble Averages

## Overview

In Section 2.1, we introduced the idea that "macroscopic properties are expected values of microscopic properties," where expected values account for the probabilities of different microstates. In Section 2.2, we derived the probability (density function) of finding a closed system in a particular microstate. In this section, we will combine these concepts to calculate macroscopic properties from expected values of microscopic properties, which we will refer to as **ensemble averages**.

## Thermodynamic Equilibrium

In Section 1.1, we defined **thermodynamic equilibrium** as "a state of simultaneous mechanical, thermal, and chemical equilibrium." We will now define mechanical, thermal, and chemical equilibrium.

```{glossary}
Thermal contact
    A state in which two systems in contact can exchange energy.

Chemical contact
    A state in which two systems in contact can exchange matter.

Mechanical equilibrium
    A state in which the net force on all particles in a system is zero.

Thermal equilibrium
    A state in which the net exchange of energy between all systems in thermal contact is zero.

Chemical equilibrium
    A state in which the net exchange of matter between all systems in chemical contact is zero.
```

## Internal Energy

The internal energy $U$ of a macroscopic system in thermodynamic equilibrium is defined as the ensemble average $\langle \cdot \rangle$ of the total energy $E$ of the system.

```{math}
U \overset{\text{def}}{=} \langle E \rangle = \sum_{i = 1}^M E_i p_i
```

where $E_i$ is the energy of microstate $i$ and $p_i$ is the probability of finding the system in microstate $i$. Substituting the probability (density function) for the canonical ensemble (i.e., the canonical distribution) from Section 2.2, we get

```{math}
:label: eq:ensemble-average-internal-energy-sum
U = \sum_{i = 1}^M E_i \frac{e^{-\beta E_i}}{Q} = \frac{1}{Q} \sum_{i = 1}^M E_i e^{-\beta E_i}
```

where $Q$ is the partition function for the canonical ensemble (i.e., the canonical partition function) from Section 2.2.

```{admonition} Why can we factor out the partition function?
:class: dropdown
Since the partition function is defined at constant number of particles, volume, and temperature, it does not depend on the (microstate) index of summation. Therefore, the partition function is a constant and can be factored out of the sum.
```

### Partial Derivative of the Partition Function with Respect to $\beta$

If we take the partial derivative of the partition function with respect to $\beta$, we get an interesting result.

```{math}
\left( \frac{\partial Q}{\partial \beta} \right)_{N, V} = \left( \frac{\partial}{\partial \beta} \sum_{i = 1}^M e^{-\beta E_i} \right)_{N, V} = \sum_{i = 1}^M \left( \frac{\partial}{\partial \beta} e^{-\beta E_i} \right)_{N, V} = -\sum_{i = 1}^M E_i e^{-\beta E_i}
```

Substituting this result into Equation {eq}`eq:ensemble-average-internal-energy-sum`, we get

```{math}
:label: eq:ensemble-average-internal-energy
U = - \frac{1}{Q} \left( \frac{\partial Q}{\partial \beta} \right)_{N, V} = -\left( \frac{\partial \ln Q}{\partial \beta} \right)_{N, V}
```

````{admonition} Where did the natural logarithm come from?
:class: dropdown
Applying the chain rule to the partial derivative of the natural logarithm of the partition function with respect to $\beta$ gives

```{math}
\left( \frac{\partial \ln Q}{\partial \beta} \right)_{N, V} = \left( \frac{\partial \ln Q}{\partial Q} \frac{\partial Q}{\partial \beta} \right)_{N, V} = \frac{1}{Q} \left( \frac{\partial Q}{\partial \beta} \right)_{N, V}
```
````

### Internal Energy of a Two-State System

In Section 2.2, we derived the canonical partition function for a two-state system as

```{math}
Q_{\text{two-state}} = e^{-\beta E_1} + e^{-\beta E_2}
```

Substituting this result into Equation {eq}`eq:ensemble-average-internal-energy`, we get

```{math}
U_{\text{two-state}} = -\left( \frac{\partial \ln Q_{\text{two-state}}}{\partial \beta} \right)_{N, V} = -\left[ \frac{\partial}{\partial \beta} \ln \left( e^{-\beta E_1} + e^{-\beta E_2} \right) \right]_{N, V}
```

Applying the chain rule to the partial derivative, we get

```{math}
U_{\text{two-state}} = \frac{E_1 e^{-\beta E_1} + E_2 e^{-\beta E_2}}{e^{-\beta E_1} + e^{-\beta E_2}}
```

````{admonition} Consistency Check
:class: dropdown
Since the denominator is $Q_{\text{two-state}}$, we can write $U_{\text{two-state}}$ as

```{math}
U_{\text{two-state}} = \frac{E_1 e^{-\beta E_1} + E_2 e^{-\beta E_2}}{Q_{\text{two-state}}}
```

Applying the definition of the canonical distribution, we get

```{math}
U_{\text{two-state}} = E_1 p_1 + E_2 p_2 = \sum_{i = 1}^2 E_i p_i = \langle E \rangle,
```

which demonstrates that Equation {eq}`eq:ensemble-average-internal-energy` is consistent with the definition of an ensemble average.
````

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, eV
from scipy.differentiate import derivative
from labellines import labelLines
from myst_nb import glue
from matplotlib.patches import Rectangle

k_B = k / eV  # Boltzmann constant in eV/K

# Define the partition function for a two-state system
def partition_function_two_state(E1, E2, T):
    beta = 1 / (k_B * T)
    return np.exp(-beta * E1) + np.exp(-beta * E2)

# Calculate the partition function for a two-state system
E1 = 0
E2 = 0.01  # Energy difference between the two states in eV
T_values = np.linspace(1, 1000, 1000)
Q_values = [partition_function_two_state(E1, E2, T) for T in T_values]

# Calculate the internal energy for a two-state system
beta_values = 1 / (k_B * T_values)
ln_Q_values = np.log(Q_values)
U_values = -np.gradient(ln_Q_values, beta_values)

# Plot the internal energy as a function of temperature
fig, ax = plt.subplots(figsize=(4, 4))
ax.plot(T_values, U_values, 'k-')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('$U_{\\text{two-state}} - E_1$ (eV)')
ax.grid(True)
ax.annotate(
    '$U_{\\text{two-state}} \\rightarrow E_1$', xy=(40, U_values[0] + 0.0001), xytext=(300, 0.001),
    arrowprops=dict(arrowstyle='->', color='b'),
    bbox=dict(boxstyle='round,pad=0.3', fc='w', ec='b'),
    ha='center', va='center', color='b'
)
x_values_high_T = np.linspace(0, 1000, 1001)
y_values_high_T = ((E1 + E2) / 2) * np.ones_like(x_values_high_T)
line = ax.plot(x_values_high_T, y_values_high_T, 'm--', label='$U_{\\text{two-state}} \\rightarrow (E_1 + E_2) / 2$')
labelLines(line, zorder=2.5)
ax.set_xlim(0, 1000)
ax.set_ylim(0, 0.01)

plt.tight_layout()
glue("fig:ensemble-average-internal-energy-two-state", fig, display=False)
plt.close(fig)
```

```{glue:figure} fig:ensemble-average-internal-energy-two-state
:figwidth: 100%
:align: center

Internal energy of a two-state system as a function of temperature. The energy difference between the two states is 0.01 eV.
```

## Heat Capacity at Constant Volume

The heat capacity at constant volume $C_V$ (to be discussed in Module 3) is defined as the partial derivative of the internal energy $U$ with respect to temperature $T$ at constant number of particles $N$ and volume $V$.

```{admonition} What Does Heat Capacity Describe?
:class: dropdown

In Section 1.1, heat was defined as "energy transferred because of a temperature difference." Heat capacity describes the amount of heat required to change the temperature of a system.
```

```{math}
C_V \overset{\text{def}}{=} \left( \frac{\partial U}{\partial T} \right)_{N, V}
```

Substituting Equation {eq}`eq:ensemble-average-internal-energy`, we get

```{math}
C_V = -\left( \frac{\partial^2 \ln Q}{\partial \beta^2} \right)_{N, V}
```
