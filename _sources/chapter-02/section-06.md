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

# 2.6. Harmonic Oscillator

## Overview

This section discusses the quantum mechanical treatment of a harmonic oscillator, along with the derivation of key thermodynamic properties using the partition function. We begin with a review of the one-dimensional case, generalize to three dimensions, and then show how the classical partition function emerges in the high-temperature limit.

## Review of the One-Dimensional Harmonic Oscillator

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.constants import k, eV
from scipy.special import hermite, factorial
from labellines import labelLines
from myst_nb import glue
from numpy.linalg import eigh

fig, axs = plt.subplot_mosaic([[0]], figsize=(4, 4))

def harmonic_potential(x, m, omega):
    """
    Calculates the potential matrix for a particle in a one-dimensional harmonic potential
    
    Parameters
    ----------
    x : array-like
        The positions of the particles.
    m : float
        The mass of the particle.
    omega : float
        The frequency of the harmonic potential
    
    Returns
    -------
    array
        The potential matrix.
    """
    V = 0.5 * m * omega**2 * x**2
    V_matrix = np.diag(V)
    return V_matrix

def off_diagonal_identity_matrix(n_points):
    matrix = np.zeros((n_points, n_points))
    for i in range(n_points):
        if i == 0:
            matrix[i + 1, i] = 1
        elif i == n_points - 1:
            matrix[i - 1, i] = 1
        else:
            matrix[i - 1, i] = 1
            matrix[i + 1, i] = 1
    return matrix

def laplacian_matrix(x):
    """
    Calculates the Laplacian matrix for a particle in a one-dimensional potential.
    
    Parameters
    ----------
    x : array-like
        The positions of the particles.
    
    Returns
    -------
    array
        The Laplacian matrix.
    """
    Delta_x = x[1] - x[0]
    n_points = len(x)
    off_diag = off_diagonal_identity_matrix(n_points)
    Laplacian = (1 / (Delta_x**2)) * (-2 * np.eye(n_points) + off_diag)
    return Laplacian

# Define constants in atomic units
hbar = 1
m = 1

# Harmonic potential
omega = 1

# Discretize the system
n_points = 2000
L = 40
x = np.linspace(-L/2, L/2, n_points)

# Construct the Hamiltonian matrix
H_harm = -hbar**2 / (2 * m) * laplacian_matrix(x) + harmonic_potential(x, m=1, omega=1)

# Solve for eigenvalues and eigenfunctions
E_harm, psi_harm = np.linalg.eigh(H_harm)

# Plot a one-dimensional quadratic potential
axs[0].plot(x, np.diag(harmonic_potential(x, m, omega)), "k-", zorder=2.5, lw=4)

# Plot the energy levels (blue lines)
for n in range(4):
    energy = axs[0].plot(x, np.ones_like(x) * E_harm[n], color='blue', label=r'$E_{%d}$' % n)
    labelLines(energy, xvals=[-L/24], zorder=2.5)

# Plot the wavefunctions (red curves)
for n in range(4):
    wavefunction = axs[0].plot(x, psi_harm[:, n] * 4 + E_harm[n], color='red', label=r'$\psi_{%d}$' % n)
    labelLines(wavefunction, xvals=[L/24], zorder=2.5)

# Format plot
axs[0].set_xlim(-L/8, L/8)
ymin = np.min(np.diag(harmonic_potential(x, m, omega)))
ymax = np.max(psi_harm[:, -1]) + E_harm[3]
yrange = ymax - ymin
ybuffer = 0.05 * yrange
axs[0].set_ylim(ymin - ybuffer, ymax + ybuffer * 2)
axs[0].set_xlabel("Position")
axs[0].set_ylabel("Energy")
axs[0].set_xticks([])
axs[0].set_yticks([])
axs[0].spines['top'].set_visible(False)
axs[0].spines['right'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)

plt.tight_layout()
glue("harmonic_oscillator", fig, display=False)
plt.close(fig)
```

```{glue:figure} harmonic_oscillator
:figwidth: 100%
:align: center
A one-dimensional harmonic oscillator. The energy levels $E_n$ are shown as horizontal blue lines, while the wavefunctions $\psi_n$ are shown as red curves.
```
