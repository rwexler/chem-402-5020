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

# Kinetic Theory

## Overview

Kinetic theory explains the macroscopic properties of gases—such as pressure and temperature—by describing the microscopic motion of particles. This lecture examines the foundational assumptions, derivation of fundamental equations, and the physical significance of kinetic theory.

---

## Foundational Assumptions of Kinetic Theory

1. **Large Number of Particles**:
   A gas consists of a statistically large number of identical particles in isotropic random motion.

2. **Point Particles**:
   The size of the particles is negligible compared to the average distance between them.

3. **Elastic Collisions**:
   Collisions between particles and with the container walls (a type of boundary) are perfectly elastic, conserving both momentum and kinetic energy.

4. **No Interparticle Forces**:
   Particles do not exert forces on one another, except during collisions.

5. **Classical Mechanics Applies**:
   Particle motion follows Newton's second law:

   ```{math}
   \vec{F} = \frac{d\vec{p}}{dt},
   ```

   where $\vec{F}$ is the net force acting on a particle, $\vec{p}$ is the linear momentum of the particle, and $t$ is time.

---

## Deriving Pressure from Particle-Wall Collisions

```{code-cell} ipython3
:tags: [hide-input]

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from myst_nb import glue

def plot_container_2d(offset=0.2):
    """Plot a 2D version of the gas particle in a container."""
    fig, ax = plt.subplots(figsize=(12, 4))

    # Draw the container as a rectangle
    Lx, Lz = 10, 2
    ax.plot([0, Lx, Lx, 0, 0], [0, 0, Lz, Lz, 0], color='black')
    ax.fill_between([0, Lx], 0, Lz, color='lightgray')

    # Draw the gas particle as a circle
    ax.plot(0.5 * Lx, 0.75 * Lz, 'o', color='blue', markersize=20, zorder=10)
    ax.text(0.5 * Lx, 0.75 * Lz, "$m$", color='white', ha='center', va='center', zorder=20, fontsize=12)

    # Add velocity arrows
    ax.annotate("", xy=(0.5 * Lx, 0.75 * Lz), xytext=(Lx, 0.75 * Lz), arrowprops=dict(arrowstyle="<-", color='red'))
    ax.text(Lx * 2 / 3, 0.75 * Lz + offset, "$v_x$", color='red', fontsize=12, ha='center', va='center')
    ax.annotate("", xy=(Lx, 0.5 * Lz), xytext=(0, 0.5 * Lz), arrowprops=dict(arrowstyle="<-", color='red'))
    ax.text(Lx * 5 / 6, 0.5 * Lz + offset, "$-v_x$", color='red', fontsize=12, ha='center', va='center')
    ax.annotate("", xy=(0, 0.25 * Lz), xytext=(0.5 * Lx, 0.25 * Lz), arrowprops=dict(arrowstyle="<-", color='red'))
    ax.text(0.25 * Lx, 0.25 * Lz + offset, "$v_z$", color='red', fontsize=12, ha='center', va='center')

    # Add length indicators
    ax.annotate("", xy=(0, -offset), xytext=(Lx, -offset), arrowprops=dict(arrowstyle="<->", color='black'))
    ax.text(Lx / 2, -2 * offset, "$L_x$", color='black', fontsize=12, ha='center', va='center')
    ax.annotate("", xy=(-offset, 0), xytext=(-offset, Lz), arrowprops=dict(arrowstyle="<->", color='black'))
    ax.text(-2 * offset, Lz / 2, "$L_z$", color='black', fontsize=12, ha='center', va='center')

    # Configure plot appearance
    ax.set_xlim(-1, 11)
    ax.set_ylim(-1, 3)
    ax.axis('off')

    return fig

# Call the function to plot
fig = plot_container_2d()

# Save and display
glue('gas-in-container', fig, display=False)
plt.close(fig)
```

```{glue:figure} gas-in-container
:name: fig-gas-in-container
:figwidth: 100%
:align: center

2D view of a gas particle (blue) in a 3D cuboid container (light gray), showing velocity $v_x$ (red arrows) and container lengths $L_x$ and $L_z$ (black bidirectional arrows). The length $L_y$ is not depicted, as it is perpendicular to the 2D view.
```

### Microscopic Picture of Pressure

Pressure is defined as the force exerted per unit area. At the microscopic scale, pressure arises from particles colliding with the container walls. Each collision exerts a force, and the cumulative effect of countless microscopic collisions produces macroscopic pressure.

#### Particle Momentum Change

Consider a particle of mass $m$ moving with velocity $v_x$ in the $x$-direction. When it collides elastically with a container wall, the momentum change $\Delta p_x$ is:

```{math}
\Delta p_x = 2mv_x.
```

#### Time Between Collisions

The time $\Delta t$ between consecutive collisions with the same wall is:

```{math}
\Delta t = \frac{2L_x}{v_x},
```

where $L_x$ is the length of the container in the $x$-direction.

#### Force Exerted on the Wall

````{margin}
```{note}
Force is exerted only during collisions.
```
````

The time-averaged force $F_{\text{p},x}$ exerted by one particle in the $x$-direction is:

```{math}
F_{\text{p},x} = \frac{\Delta p_x}{\Delta t} = \frac{mv_x^2}{L_x},
```

where $\text{p}$ denotes one particle.

### Total Pressure

For $N$ identical particles moving isotropically and randomly in a container of volume $V$, the pressure $P$ is:

```{math}
P = \frac{1}{V} \sum_{i=1}^N \frac{1}{3} m v_i^2,
```

where $v_i$ is the speed of the $i^{\text{th}}$ particle. Substituting the sum with an average, we obtain:

<!-- ````{important} -->
```{math}
:label: pressure-kinetic-theory
P = \frac{N m \langle v^2 \rangle}{3V}.
```

Equation {eq}`pressure-kinetic-theory` relates the macroscopic pressure to the microscopic average squared speed $\langle v^2 \rangle$ of the gas particles, which is proportional to their average kinetic energy.

---

## Kinetic Energy and Temperature

By definition, the average translational kinetic energy $\langle E_\text{kin} \rangle$ is:

```{math}
\langle E_\text{kin} \rangle = \frac{1}{2} m \langle v^2 \rangle.
```

Using the equation $PV = Nk_\text{B}T$, referred to as the ideal gas law (to be discussed in Lecture 3), we find:
<!-- Reference Lecture 3 -->

```{math}
:label: equipartition-theorem
\frac{1}{2} m \langle v^2 \rangle = \frac{3}{2} k_\text{B} T,
```

where $k_\text{B}$ is the [Boltzmann constant](https://physics.nist.gov/cgi-bin/cuu/Value?k) (i.e., $k_\text{B} = 1.381 \times 10^{-23} \, \text{J/K}$) and $T$ is the absolute temperature of the gas. Equation {eq}`equipartition-theorem`, termed the equipartition theorem, shows that the macroscopic temperature is directly proportional to the microscopic average kinetic energy of the gas particles, indicating the equal partition of energy among degrees of freedom.
<!-- Clarify degrees of freedom. -->

---

## Estimating Particle Speeds

Solving Equation {eq}`equipartition-theorem` for $\langle v^2 \rangle$ and taking the square root, we calculate the root-mean-square speed $v_{\text{rms}}$ of the gas particles:

```{math}
:label: rms-speed
v_{\text{rms}} = \sqrt{\frac{3k_\text{B}T}{m}},
```

which is a statistical measure of the average speed of the gas particles. This demonstrates that hotter gases have faster-moving particles than colder gases and that lighter particles move faster than heavier ones at the same temperature.

````{admonition} Example: Estimating the Speed of Air Molecules at Room Temperature
:class: dropdown

As N<sub>2</sub> molecules constitute 78.084% of air molecules, we can estimate the speed of air molecules at room temperature as the root-mean-square speed of N<sub>2</sub> molecules at 300 K. Using the relative molecular mass $M_\text{r}$ of N<sub>2</sub> (i.e., $M_\text{r} = 28.0134 \, \text{g/mol}$[^1]), we find:

```{math}
v_{\text{rms}} = \sqrt{\frac{3k_\text{B}T}{m}} = \sqrt{\frac{3k_\text{B}T}{M_\text{r} / N_\text{A}}} \approx 517 \, \text{m/s} \approx 1{,}156 \, \text{mph}.
```

This calculation shows that air molecules move several hundreds of miles per hour faster than a Boeing 747-8 aircraft at cruising speed (about 660 mph).[^2]

[^1]: [https://webbook.nist.gov/cgi/cbook.cgi?Formula=N2&NoIon=on&Units=SI](https://webbook.nist.gov/cgi/cbook.cgi?Formula=N2&NoIon=on&Units=SI)
[^2]: [https://www.boeing.com/commercial/747-8/design-highlights#technologically-advanced](https://www.boeing.com/commercial/747-8/design-highlights#technologically-advanced)
````

<!-- Just note that N<sub>2</sub> is actually *diatomic*. Strictly speaking, the factor of 3 in the RMS formula comes from *translational* degrees of freedom only; it remains valid for *translational speed*. Diatomic gases do have rotational (and possibly vibrational) modes, but that mainly affects *internal* energy—translation is still 3D motion. -->
