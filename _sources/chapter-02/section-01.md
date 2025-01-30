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

# Canonical Ensemble

## Overview

```{mermaid}
---
config: "{
  \"theme\": \"neutral\",
  \"look\": \"handDrawn\",
  \"layout\": \"elk\"
}"
---
flowchart LR
  %% Beginnings and ending
  CM([Classical mechanics])
  QM([Quantum mechanics])
  Th([Thermodynamics])

  %% Processes
  KT[[Kinetic theory]]
  SM[[Statistical mechanics]]

  %% Inputs/Outputs
  StateCM[/"<i>r<sup>N</sup></i>, <i>p<sup>N</sup></i>"/]
  StateQM[/"&Psi;(<i>r<sup>N</sup></i>)"/]
  
  %% Decision
  LimitCM{"Classical limit?"}

  subgraph Microscopic World
    QM
    CM
    StateCM
    StateQM
    LimitCM
  end
  
  subgraph Bridges
    KT
    SM
  end
  
  subgraph Macroscopic World
    Th
  end
  
  CM --> StateCM
  StateCM --> KT --> Th
  StateCM --> SM --> Th
  
  QM --> LimitCM
  LimitCM -- Yes --> CM
  LimitCM -- No --> StateQM --> SM
```

In chemistry, our main goal is often to **predict** and **control** the future of chemical systems. For example:

- If we mix two reagents, will they react or remain inert?
- Will they be miscible (like ethanol in water) or immiscible (like oil in water)?

**Thermodynamics** gives us powerful laws for understanding energy changes and equilibria, but these laws describe the system at a **macroscopic** level (e.g., temperature, pressure, volume, composition). **Statistical mechanics** helps us connect these macroscopic observables to the microscopic world of molecules, atoms, and their interactions. It uses the results of **quantum mechanics** (and sometimes **classical mechanics**) to interpret and predict the behavior of a large number of particles.

<!-- We will organize our discussion using five main sections ("chalkboards"), each of which addresses a crucial piece of the puzzle. -->

---

## Motivations — Connecting Microscopic States to Macroscopic Observables

### Predicting the Future of Chemical Systems

Thermodynamics tells us that the equilibrium state of a system is determined by minimizing a free energy (e.g., the Gibbs free energy at constant temperature and pressure). However, it does not directly tell us **how** molecules arrange themselves or **why** certain processes occur at the molecular level. This is where **statistical mechanics** excels.

1. **The Hard Way**: Simulate each particle quantum mechanically or classically.
   - **Quantum Dynamics**: Solve the many-body time-dependent Schrödinger equation for each electron and nucleus. This is quickly intractable for large systems.
   - **Classical Dynamics**: Solve Newton’s equations of motion for every particle in the system. This can handle larger systems (up to millions or billions of particles), but still falls short of the $\sim 10^{23}$ particles encountered in everyday laboratory-scale chemistry.

2. **The Easy Way**: Use **statistical mechanics** and **thermodynamics**:
   - Thermodynamics provides general relationships (e.g., equilibrium occurs at minimum free energy).
   - Statistical mechanics **computes** free energy from the **microscopic** properties of molecules, thereby revealing *why* the thermodynamics "rules" hold.

In essence, **statistical mechanics** bridges the gap: it starts from **microscopic states** (positions, momenta, wavefunctions, etc.) and delivers **macroscopic observables** (pressure, temperature, free energy, etc.).

---

## Probability Distributions & the Fundamental Postulate

### Probability Distributions in Statistical Mechanics

Instead of tracking exactly which microstate the system is in at every moment, we use **ensembles**—conceptual collections of many virtual copies of our system under the same macroscopic conditions. We then describe the likelihood (probability) of being in each microstate.

One of the most important assumptions in statistical mechanics is the **Fundamental "Equal A Priori Probability" Postulate**, often stated as:

> For an **isolated system** with an **exactly known energy** and **exactly known composition**, the system can be found with **equal probability** in any microstate consistent with that knowledge.

### Names & Interpretations

- **Principle of Indifference**: In the absence of any additional information, all states are equally probable.
- **Sherlock Holmes Principle**: "When you have eliminated the impossible, whatever remains, however improbable, must be the truth."

### The Ergodic Hypothesis

Closely related is the **ergodic hypothesis**:

$$
\text{time average} = \text{ensemble average}.
$$

It states that measuring an observable (e.g., energy) of a **single** system over a long time yields the same average as measuring an **ensemble** of many identical systems once. This hypothesis is generally (but not always) assumed in chemistry.

**Non-ergodic systems** (e.g., glass slowly turning into quartz over geological timescales) violate this assumption, as they do not sample all energetically accessible states on a human timescale.

---

## Microstates, Macrostates, & Degeneracy

### Microstates

A **microstate** is any unique (microscopic) configuration of a system consistent with its macroscopic constraints (e.g., total energy, number of particles, volume). Microstates can be described by:

- Classical positions and momenta of all particles: $\{(\mathbf{r}_i, \mathbf{p}_i)\}$.
- Quantum mechanical wavefunction specifying all electrons and nuclei.
- Spin configurations in magnetic systems.
- Complete quantum numbers specifying each particle’s state.

### Macrostates

A **macrostate** is specified by the **minimum set of thermodynamic variables** (e.g., $N, P, V, T$) that describe the system on a bulk scale. Many different microstates can correspond to the **same** macrostate.

Examples:

- A gas in a container described by $(N, V, T)$.
- A solution described by $(T, P, \text{composition})$.

### Degeneracy

An **energy level** may correspond to more than one distinct microstate—this number is called the **degeneracy** $\Omega$. Formally:

> $\Omega$ = number of microstates all having the same energy $E$.

This distinction is crucial:

- **Microstate**: One **specific** arrangement (or wavefunction) of the system.
- **Energy level**: A **value** of energy the system can have.
- **Degeneracy**: How many microstates share that energy.

---

## Microcanonical Ensemble & Isolated Systems

### The Microcanonical Ensemble

An **isolated** system has a boundary through which **no** matter or energy is exchanged. Such a system has exactly known:

1. **Energy** $E$
2. **Number of particles** $N$
3. **Volume** $V$ (or area $A$ in 2D)

The collection of all microstates consistent with these constraints is called the **microcanonical ensemble**. According to the fundamental postulate, each microstate in this ensemble is **equally likely**.

Mathematically, if there are $\Omega$ microstates all sharing the same $E, N, V$, then

$$
P_i = \frac{1}{\Omega}, \quad i = 1, 2, \dots, \Omega.
$$

### Example: 2D Box with Five Particles

Imagine a system of five non-interacting particles confined to a 2D box. The constraints are:

- $N=5$
- $E=E_0$
- Area $A=L^2$

Each microstate corresponds to a distinct arrangement of the five particles (and their momenta). Since the system is isolated, all these microstates must share the same total energy $E_0$. By the postulate, each is **equally** probable.

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
from myst_nb import glue

# For reproducibility
rng = np.random.default_rng(seed=42)

# Define the number of particles and the box dimension
N = 5
L = 10

# Create a 3x3 grid of subplots
fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(12, 12))

for i in range(3):
    for j in range(3):
        # Random positions in [0, L)
        x = rng.random(N) * L
        y = rng.random(N) * L

        # Random momenta centered around 0
        px = (rng.random(N) - 0.5) * 2
        py = (rng.random(N) - 0.5) * 2

        # Plot positions
        axs[i, j].scatter(x, y, c='blue', s=50, label='Positions')
        # Plot momenta as arrows
        axs[i, j].quiver(
            x, y, px, py,
            color='red',
            angles='xy',
            scale_units='xy',
            scale=1,
            alpha=0.8,
            label='Momenta'
        )

        # Set axis limits and aspect ratio
        axs[i, j].set_xlim(0, L)
        axs[i, j].set_ylim(0, L)
        axs[i, j].set_aspect('equal')

        # Remove tick labels for a cleaner look
        axs[i, j].set_xticks([])
        axs[i, j].set_yticks([])

        # Add a descriptive title
        microstate_index = 3 * i + j + 1
        axs[i, j].set_title(
            f"Microstate {microstate_index}\n"
            f"Macrostate: $N = {N}, E = E_0, A = L^2$"
        )

# Add a main title
plt.suptitle("Microcanonical Ensemble: 9 Example Microstates", fontsize=16)
fig.tight_layout()

# Save figure for later embedding with myst_nb's glue
glue("microcanonical_ensemble", fig, display=False)

# Close figure to prevent duplicate display
plt.close(fig)
```

```{glue:figure} microcanonical_ensemble
:name: microcanonical-ensemble
:figwidth: 100%
:align: center

**Microcanonical Ensemble**: This figure shows nine representative microstates of an isolated 2D system containing five **classical** (i.e., well-defined positions and momenta) and non-interacting particles. All nine microstates share the same macrostate constraints—namely $N = 5$, $E = E_0$, and a box area $A = L^2$—yet differ in the precise distribution of positions and momenta among the particles.
```

---

## Linking Thermodynamics & Statistics

### A Simple Example: Average Number of Particles

In an isolated system, $N$ is constant. Therefore, the **average number of particles** $\langle N\rangle$ trivially equals $N$. Formally:

$$
\langle N \rangle
= \sum_{i=1}^{\Omega} N_i \, P_i
= \sum_{i=1}^{\Omega} 5 \cdot \frac{1}{\Omega}
= 5.
$$

Although this result is obvious, it demonstrates the basic approach in statistical mechanics: a **weighted average** over all microstates, with probabilities given by $P_i$.

### Other Averages

By the same principle, we can compute **average energy** $\langle E\rangle$ or **average volume** $\langle V\rangle$. For an isolated system, these remain constant, so the average equals the known value. In more complex systems (e.g., a **canonical ensemble** exchanging heat but not matter with surroundings), the probability distribution is no longer uniform, and $\langle E\rangle$ can fluctuate in more interesting ways.

### The Big Picture

1. **Statistical Mechanics**:
   - Tells us how to compute probabilities of different microstates.
   - Links microscopic details to macroscopic observables (e.g., $T, P, \langle E\rangle$).
2. **Thermodynamics**:
   - Summarizes the macroscopic laws (e.g., equilibrium at minimum free energy).
   - Provides relationships such as $dU = TdS - PdV + \dots$.
3. **Quantum/ Classical Mechanics**:
   - Determines the energy levels or potential energy surfaces relevant to each microstate.

Putting them all together allows us to predict **equilibrium constants**, **reaction rates** (with the help of kinetics), **phase behavior** (miscibility, vapor pressure, etc.), and **countless** other chemical phenomena.

---

## Summary

1. **Motivations**: We need to connect the *microscopic* (positions, momenta, wavefunctions) to the *macroscopic* (energy, pressure, volume, composition) to predict and understand chemical phenomena.
2. **Probability Distributions & Fundamental Postulate**: In an isolated system of fixed $E$ and $N$, each accessible microstate is equally probable. This is the cornerstone of statistical mechanics.
3. **Microstates, Macrostates, & Degeneracy**: A single macrostate corresponds to many possible microstates. Degeneracy $\Omega$ counts how many microstates share a given energy.
4. **Microcanonical Ensemble & Isolated Systems**: All microstates with the same $E, N, V$ are equally likely. This is the simplest ensemble in statistical mechanics.
5. **Linking Thermodynamics & Statistics**: Averages of thermodynamic quantities emerge from probability-weighted sums over microstates. This is how macroscopic observables arise from microscopic theory.
