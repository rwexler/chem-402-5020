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

# 2.4. Molecular Partition Functions

## Overview

In Sections 2.2 and 2.3, we derived the canonical partition function and the ensemble averages for systems of one particle. In this section, we extend these concepts to closed systems of **many identical and independent particles with internal degrees of freedom**. We will also introduce the concept of **indistinguishability**, which is central to the statistical mechanics of quantum systems.

```{note}
**Assumption (Single Occupancy):**
Throughout this section, we assume that each one‐particle microstate can be occupied by **at most one** particle. Under this assumption, the total number of ways to place $N$ distinguishable particles in $M$ one‐particle microstates is $\tfrac{M!}{(M-N)!}$ rather than $M^N$. If multiple occupancy were allowed, the counting would differ accordingly (e.g., $M^N$ for unlimited occupancy).
```

## Partition Function for Distinguishable Particles

````{admonition} Two Distinguishable Particles in a Four-State System
:class: dropdown

The $N$-particle microstates of a system of two distinguishable particles in a four-state system are given by the following table:

```{list-table}
:header-rows: 1

* - N-Particle Microstate
  - 1-Particle Microstate 1
  - 1-Particle Microstate 2
  - 1-Particle Microstate 3
  - 1-Particle Microstate 4
* - I
  - a
  - b
  - 
  - 
* - II
  - a
  - 
  - b
  - 
* - III
  - a
  - 
  - 
  - b
* - IV
  - 
  - a
  - b
  - 
* - V
  - 
  - a
  - 
  - b
* - VI
  - 
  - 
  - a
  - b
* - VII
  - b
  - a
  - 
  - 
* - VIII
  - b
  - 
  - a
  - 
* - IX
  - b
  - 
  - 
  - a
* - X
  - 
  - b
  - a
  - 
* - XI
  - 
  - b
  - 
  - a
* - XII
  - 
  - 
  - b
  - a
```

\
The number of $N$-particle microstates of this system is $12$.

Here, we treat each of the four one‐particle states as an exclusive "slot," so two distinguishable particles must occupy different states, yielding $4 \times 3 = 12$ arrangements rather than $4 \times 4 = 16$.

---

In general, the number of $N$-particle microstates of a system of $N$ distinguishable particles in a system with $M$ one-particle microstates is given by the number of permutations of $N$ particles in $M$ states:

```{math}
\frac{M!}{(M - N)!}
```

This assumes no more than one particle per one‐particle state.
````

The total energy $E$ of a system of $N$ distinguishable and independent particles is given by the sum of the energies $\varepsilon$ of each particle ($a$, $b$, $c$, etc.):

```{math}
E_l = \varepsilon_i^a + \varepsilon_j^b + \varepsilon_k^c + \cdots
```

where $i$, $j$, and $k$ are the indices of the (*one-particle*) microstates of the particles $a$, $b$, and $c$, respectively, and $l$ is the index of the (*many-particle* or *N-particle*) microstate of the system.

The canonical partition function $Q$ for a system of $N$ distinguishable and independent particles is given by the sum over all possible microstates:

```{math}
Q = \sum_l e^{-\beta E_l} = \sum_{i = 1}^M \sum_{j = 1}^M \sum_{k = 1}^M \cdots e^{-\beta (\varepsilon_i^a + \varepsilon_j^b + \varepsilon_k^c + \cdots)} = \sum_{i = 1}^M e^{-\beta \varepsilon_i^a} \sum_{j = 1}^M e^{-\beta \varepsilon_j^b} \sum_{k = 1}^M e^{-\beta \varepsilon_k^c} \cdots = \boxed{q_a^N q_b^N q_c^N \cdots}
```

where $q_a$, $q_b$, and $q_c$ are the canonical partition functions for the particles $a$, $b$, and $c$, respectively.

If the particles are identical, we can write the canonical partition function as

```{math}
Q = q^N
```

where $q$ is the canonical partition function for one particle.

## Partition Function for Indistinguishable Particles

````{admonition} Two Indistinguishable Particles in a Four-State System
:class: dropdown

The $N$-particle microstates of a system of two indistinguishable particles in a four-state system are given by the following table:

```{list-table}
:header-rows: 1

* - N-Particle Microstate
  - 1-Particle Microstate 1
  - 1-Particle Microstate 2
  - 1-Particle Microstate 3
  - 1-Particle Microstate 4
* - I
  - x
  - x
  - 
  - 
* - II
  - x
  - 
  - x
  - 
* - III
  - x
  - 
  - 
  - x
* - IV
  - 
  - x
  - x
  - 
* - V
  - 
  - x
  - 
  - x
* - VI
  - 
  - 
  - x
  - x
```

\
The number of $N$-particle microstates of this system is $6$.

---

In general, the number of $N$-particle microstates of a system of $N$ indistinguishable particles in a system with $M$ one-particle microstates is given by the number of combinations of $N$ particles in $M$ states:

```{math}
\frac{1}{N!} \frac{M!}{(M - N)!}
```

This assumes no more than one particle per one‐particle state.
````

If the particles are indistinguishable and identical, we can write the canonical partition function as

```{math}
Q = \frac{q^N}{N!}
```

## Partition Function for Molecules

Within the Born-Oppenheimer approximation, the total energy of a molecule is given by the sum of the translational, rotational, vibrational, and electronic energies:

```{math}
\varepsilon_{\lambda} = \varepsilon_i^{\text{trans}} + \varepsilon_j^{\text{rot}} + \varepsilon_k^{\text{vib}} + \varepsilon_l^{\text{elec}}
```

where $i$, $j$, $k$, and $l$ are the indices of the (*one-degree-of-freedom*) microstates of the translational, rotational, vibrational, and electronic energies, respectively, and $\lambda$ is the index of the (*many-degree-of-freedom*) microstate of the molecule.

The canonical partition function $q$ for a molecule is given by the sum over all possible microstates:

```{math}
\begin{aligned}
q &= \sum_{\lambda} e^{-\beta \varepsilon_{\lambda}} = \sum_{i} \sum_{j} \sum_{k} \sum_{l} e^{-\beta (\varepsilon_i^{\text{trans}} + \varepsilon_j^{\text{rot}} + \varepsilon_k^{\text{vib}} + \varepsilon_l^{\text{elec}})} = \sum_{i} e^{-\beta \varepsilon_i^{\text{trans}}} \sum_{j} e^{-\beta \varepsilon_j^{\text{rot}}} \sum_{k} e^{-\beta \varepsilon_k^{\text{vib}}} \sum_{l} e^{-\beta \varepsilon_l^{\text{elec}}} \\
&= q_{\text{trans}} q_{\text{rot}} q_{\text{vib}} q_{\text{elec}}
\end{aligned}
```

where $q_{\text{trans}}$, $q_{\text{rot}}$, $q_{\text{vib}}$, and $q_{\text{elec}}$ are the canonical partition functions for the translational, rotational, vibrational, and electronic energies, respectively.
