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

# 2.8. Molecular Statistical Mechanics

## Overview

## Symmetry of Linear Rigid Rotors

## Nonlinear Rigid Rotors

## Vibrations of Polyatomic Molecules

## Electronic Partition Function

## Summary

```{list-table} Partition functions for various systems. Repeated entries are denoted by a down arrow.
:header-rows: 1
:name: tab:partfunc

* - System
  - $q_{\text{trans}}$
  - $q_{\text{rot}}$
  - $q_{\text{vib}}$
  - $q_{\text{elec}}$
* - Atom
  - $\dfrac{V}{\Lambda^3}$
  - 
  - 
  - $g_1$
* - Diatomic
  - $\downarrow$
  - $\dfrac{T}{\sigma\,\Theta_{\text{rot}}}$
  - $\dfrac{e^{-\Theta_{\text{vib}}/(2T)}}{1-e^{-\Theta_{\text{vib}}/T}}$
  - $\downarrow$
* - **Polyatomic**
  - 
  - 
  - 
  - 
* - *Linear*
  - $\downarrow$
  - $\downarrow$
  - $\displaystyle\prod_{j=1}^\alpha \dfrac{e^{-\Theta_{\text{vib},j}/(2T)}}{1-e^{-\Theta_{\text{vib},j}/T}}$
  - $\downarrow$
* - *Nonlinear*
  - 
  - 
  - 
  - 
* - Spherical
  - $\downarrow$
  - $\dfrac{\sqrt{\pi}}{\sigma}\left(\dfrac{T}{\Theta_{\text{rot}}}\right)^{3/2}$
  - $\downarrow$
  - $\downarrow$
* - Symmetric
  - $\downarrow$
  - $\dfrac{\sqrt{\pi}}{\sigma}\left(\dfrac{T}{\Theta_{\text{rot,A}}}\right)\left(\dfrac{T}{\Theta_{\text{rot,C}}}\right)^{1/2}$
  - $\downarrow$
  - $\downarrow$
* - Asymmetric
  - $\downarrow$
  - $\dfrac{\sqrt{\pi}}{\sigma}\left(\dfrac{T^3}{\Theta_{\text{rot,A}}\,\Theta_{\text{rot,B}}\,\Theta_{\text{rot,C}}}\right)^{1/2}$
  - $\downarrow$
  - $\downarrow$
```
