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

# 1.1. Course Introduction

## Overview

This course covers the fundamental concepts of thermodynamics and statistical mechanics, along with their applications to chemical systems. It builds upon the thermodynamic principles introduced in [Chem 106](https://chemistry.wustl.edu/chemistry-105-106) and [Chem 112](https://chemistry.wustl.edu/chemistry-111-112), and the quantum mechanics taught in [Chem 105](https://chemistry.wustl.edu/chemistry-105-106), [Chem 111](https://chemistry.wustl.edu/chemistry-111-112), and [Chem 401](https://chemistry.wustl.edu/physical-chemistry-i). By connecting molecular-level behavior to macroscopic thermodynamic observations, students will learn how theory underpins real-world chemical processes.

---

## Why Study Thermodynamics and Statistical Mechanics?

### Bridging the Microscopic and Macroscopic Worlds

Typical chemical systems contain on the order of [Avogadro's number](https://physics.nist.gov/cgi-bin/cuu/Value?na) of particles (i.e., $N_\text{A} = 6.022 \times 10^{23}$). Thermodynamics abstracts the complexity of such large collections of molecules into a framework for predicting:

- Reaction spontaneity,
- Equilibrium states,
- Phase transitions, and more.

### Modern Applications

Thermodynamics and statistical mechanics are central to diverse fields, including industrial chemistry, materials science, and biochemistry. The examples below highlight just a few applications:

::::{grid}
:gutter: 3

:::{grid-item-card} Industrial Chemistry
:text-align: center
<div class="center-text">
Predicting reaction spontaneity and modeling large-scale chemical processes.
</div>

```{figure} ../_static/section-01/industrial_chemistry.jpg
---
name: industrial-chemistry
---
Hydrogen production via steam–methane reforming.[^1]

[^1]: [https://www.energy.gov/eere/fuelcells/hydrogen-production-natural-gas-reforming](https://www.energy.gov/eere/fuelcells/hydrogen-production-natural-gas-reforming)
```

:::

:::{grid-item-card} Materials Chemistry
:text-align: center
<div class="center-text">
Designing advanced materials for optical, electronic, or mechanical applications.
</div>

```{figure} ../_static/section-01/materials_chemistry.jpg
---
name: materials-chemistry
---
Electricity generation using solar panels.[^2]

[^2]: [https://www.energy.gov/eere/solar/cadmium-telluride](https://www.energy.gov/eere/solar/cadmium-telluride)
```

:::

:::{grid-item-card} Biochemistry
:text-align: center
<div class="center-text">
Investigating protein folding and enzyme activity in biological systems.
</div>

```{figure} ../_static/section-01/biochemistry.jpg
---
name: biochemistry
---
Predicting protein structures with AlphaFold.[^3]

[^3]: [https://alphafold.ebi.ac.uk/entry/P21852](https://alphafold.ebi.ac.uk/entry/P21852)
```

:::
::::

---

## Key Definitions

### Thermodynamic Systems

```{code-cell} ipython3
:tags: [hide-input]

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from myst_nb import glue

# Helper function to plot a system
def plot_system(ax, title, annotations, boundary_color='b'):
    box = mpatches.FancyBboxPatch((0, 0), 1, 1, boxstyle='roundtooth', ec=boundary_color, fc='w')
    ax.add_patch(box)
    ax.set_title(title, fontsize=14)
    ax.text(0.5, 0.5, 'System', ha='center', va='center', fontsize=12)
    ax.text(0.5, -0.65, 'Surroundings', ha='center', va='center', fontsize=12)
    ax.text(0.5, 1.3, 'Boundary', ha='center', va='bottom', fontsize=12, color=boundary_color)
    for annotation in annotations:
        if "arrowprops" in annotation:  # Arrow annotations
            ax.annotate('', **annotation)
        else:  # Text annotations
            ax.text(**annotation)
    ax.set_xlim(-1, 2)
    ax.set_ylim(-1, 2)
    ax.set_aspect('equal')
    ax.axis('off')

# Define annotations for each system
annotations = [
    [],  # Isolated system (no arrows)
    [  # Closed system (energy arrow)
        dict(xy=(-0.6, 0.15), xytext=(0.15, 0.15), arrowprops=dict(arrowstyle='<->', color='r')),
        dict(x=-1, y=0.3, s='Energy', ha='left', va='bottom', fontsize=12, color='r'),
    ],
    [  # Open system (energy + matter arrows)
        dict(xy=(-0.6, 0.15), xytext=(0.15, 0.15), arrowprops=dict(arrowstyle='<->', color='r')),
        dict(xy=(0.85, 0.15), xytext=(1.6, 0.15), arrowprops=dict(arrowstyle='<->', color='m')),
        dict(x=-1, y=0.3, s='Energy', ha='left', va='bottom', fontsize=12, color='r'),
        dict(x=2, y=0.3, s='Matter', ha='right', va='bottom', fontsize=12, color='m'),
    ],
]

titles = ["Isolated system", "Closed system", "Open system"]

fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for i, ax in enumerate(axes):
    plot_system(ax, titles[i], annotations[i])

glue('thermo-systems', fig, display=False)
plt.close(fig)
```

```{glue:figure} thermo-systems
:name: thermo-systems
:figwidth: 100%
:align: center

**Types of thermodynamic systems:** (a) **Isolated**—no exchange of energy or matter; (b) **Closed**—exchanges energy but not matter; and (c) **Open**—exchanges both energy and matter.
```

```{glossary}
System
  The portion of the universe chosen for study, separated from its surroundings by a boundary.

Surroundings
  Everything external to the system that can exchange energy or matter with it.

Boundary
  The interface separating a system from its surroundings.

Isolated system
  Exchanges neither energy nor matter with its surroundings.

Closed system
  Exchanges energy but not matter with its surroundings.

Open system
  Exchanges both energy and matter with its surroundings.
```

### State of a System

```{mermaid}
---
config: "{
  \"flowchart\": {
    \"htmlLabels\": false
  },
  \"theme\": \"neutral\",
  \"look\": \"handDrawn\",
  \"layout\": \"elk\"
}"
---
flowchart TB
  subgraph B["Microscopic"]
    direction TB
      B1["Classical"] ~~~ B2["Quantum"]
  end

  subgraph C["Macroscopic"]
    direction TB
      C1["Equilibrium"]
  end

  subgraph C1["Equilibrium"]
    direction TB
      C11["Thermodynamic"]
  end

  subgraph C11["Thermodynamic"]
    direction LR
      C111["State variables"] ~~~ C112["State functions"] & C113["Path functions"]
  end

  subgraph C111["State variables"]
    direction TB
      C1111["Extensive"] ~~~ C1112["Intensive"]
  end

  A["State of a system"] --> B & C
```

```{glossary}
Particle
  A microscopic entity such as an atom, molecule, or ion.

Microscopic state (classical)
  Positions and momenta of all particles in the system.

Microscopic state (quantum)
  The wavefunction describing the system's particles.

Equilibrium
  A condition in which macroscopic properties remain constant over time.

Thermodynamic equilibrium
  A state of simultaneous mechanical, thermal, and chemical equilibrium.

Thermodynamic state
  A set of macroscopic variables defining a system in equilibrium.

State variable
  A property that defines a system's state.

State function
  A property depending only on the system's state, not on the path taken.

Equation of state
  A mathematical relationship among state variables.

Path function
  A property depending on the process or path taken between states.

Process
  A transformation changing a system from one state to another.

Extensive property
  A property proportional to system size (e.g., volume, entropy).

Intensive property
  A property independent of system size (e.g., temperature, pressure).
```

---

## Basic Forms of Energy and Energy Transfer

```{mermaid}
---
config: "{
  \"flowchart\": {
    \"htmlLabels\": false
  },
  \"theme\": \"neutral\",
  \"look\": \"handDrawn\",
  \"layout\": \"elk\"
}"
---
flowchart TB
  subgraph H1["Electrostatic"]
    direction TB
      H11["Charge–charge"] ~~~ H12["Charge–dipole"]
  end
  subgraph H2["van der Waals"]
    direction TB
      H21["Dipole–dipole"] ~~~ H22["Induced dipole"] ~~~ H23["Dispersion"]
  end
  subgraph H["Electric"]
    direction LR
      H1 ~~~ H2
  end
  subgraph I1["Intramolecular"]
    direction TB
      I11["Covalent"] ~~~ I12["Ionic"] ~~~ I13["Metallic"]
  end
  subgraph I2["Intermolecular"]
    direction TB
      I21["Hydrogen"] ~~~ I22["van der Waals"]
  end
  subgraph I["Chemical"]
    direction LR
      I1 ~~~ I2
  end
  subgraph J["Mechanical"]
    direction TB
      J1["Expansion"] ~~~ J2["Surface expansion"] ~~~ J3["Extension"]
  end
  subgraph G["Work"]
    direction LR
      J ~~~ G1["Electrical"]
  end

  A["Basic forms of energy and energy transfer"] --> B["Energy"] & C["Energy transfer"]
  B --> D["Kinetic"] & E["Potential"]
  C --> F["Heat"] & G
  E --> H & I
  H1 & H2
  I1 & I2

  style H11 stroke:red, stroke-width:2px
  style H2 stroke:blue, stroke-width:2px
  style I12 stroke:red, stroke-width:2px
  style I22 stroke:blue, stroke-width:2px
```

### Energy

```{glossary}
Kinetic energy
  Energy due to motion (e.g., a moving particle).

Potential energy
  Energy due to position or configuration (e.g., a stretched spring).

Internal energy
  The total microscopic kinetic and potential energy of a system, averaged over its microstates.
```

### Energy Transfer

```{glossary}
Work
  Energy transferred when a force acts over a distance (e.g., lifting a mass).

Heat
  Energy transferred because of a temperature difference (e.g., conduction from hot to cold).
```

---

## Important Units

### SI Units

```{list-table} **Base SI Units**
:header-rows: 1
:name: physical-quantities-basic

* - Quantity
  - Unit
  - Symbol
* - Time
  - second
  - s
* - Length
  - meter
  - m
* - Mass
  - kilogram
  - kg
* - Temperature
  - kelvin
  - K
```

[^4]

[^4]: [https://www.nist.gov/pml/owm/metric-si/si-units](https://www.nist.gov/pml/owm/metric-si/si-units)

```{list-table} **Derived SI Units**
:header-rows: 1
:name: physical-quantities-units

* - Quantity
  - Unit
  - Symbol
  - Conversion
* - Frequency
  - hertz
  - Hz
  - $1 \,\text{Hz} = 1 \,\text{s}^{-1}$
* - Force
  - newton
  - N
  - $1 \,\text{N} = 1 \,\text{kg m s}^{-2}$
* - Pressure
  - pascal
  - Pa
  - $1 \,\text{Pa} = 1 \,\text{N m}^{-2}$
* - Energy
  - joule
  - J
  - $1 \,\text{J} = 1 \,\text{N m}$
```

### Non-SI Units

```{list-table} **Non-SI Units**
:header-rows: 1
:name: non-si-units

* - Quantity
  - Unit
  - Symbol
  - Conversion
* - Pressure
  - bar
  - bar
  - $1 \,\text{bar} = 1 \times 10^5 \,\text{Pa}$
* - 
  - atmosphere
  - atm
  - $1 \,\text{atm} \approx 1.01325 \,\text{bar}$
* -
  - torr
  - torr
  - $1 \,\text{torr} = \frac{1}{760}\,\text{atm}$
* -
  - millimeters of mercury
  - mmHg
  - $1 \,\text{mmHg} = 1 \,\text{torr}$
* - Energy
  - electronvolt
  - eV
  - $1 \,\text{eV} = 1.602 \times 10^{-19} \,\text{J}$
* -
  - calorie
  - cal
  - $1 \,\text{cal} = 4.184 \,\text{J}$
```

[^5] [^6] [^7]

[^5]: [https://www.nist.gov/pml/sensor-science/thermodynamic-metrology/unit-conversions](https://www.nist.gov/pml/sensor-science/thermodynamic-metrology/unit-conversions)  
[^6]: [https://physics.nist.gov/cgi-bin/cuu/Value?evj](https://physics.nist.gov/cgi-bin/cuu/Value?evj)  
[^7]: [https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8#C](https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8#C)
