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

# 3.2. Applications of the First Law

## Overview

## Thermodynamic Processes

<!-- - Irreversible -->
<!-- - Reversible -->
<!-- - Quasi-static -->

```{glossary}
Quasi-static process
  A process carried out *infinitesimally slowly*, ensuring the system remains *near equilibrium* at each step. Any deviations from equilibrium are so small (infinitesimal) that each intermediate state is well-defined thermodynamically.

Reversible process
  An *idealized* process carried out *quasistatically* and with *no dissipative effects*, so that it can be undone without leaving any net change in the system or surroundings.
```

## How to Apply the First Law

1. Choose two independent variables
2. Rewrite the first law
3. Apply constraints
4. Define the system
5. Integrate the first law

## Example: $V$ and $T$ as Independent Variables

### Rewrite the First Law

```{math}
\begin{aligned}
\delta q &= dU - \delta w \\
&= dU + PdV \\
&= \left( \frac{\partial U}{\partial T} \right)_V dT + \left( \frac{\partial U}{\partial V} \right)_T dV + PdV \\
&= C_V dT + \left[ \left( \frac{\partial U}{\partial V} \right)_T + P \right] dV
\end{aligned}
```

### Apply Constraints

| Constraint | Equation | First Law |
|------------|----------|-----------|
| Isochoric  | $dV = 0$ | $\delta q = C_V dT$ |
| Isothermal | $dT = 0$ | $\delta q = \left[ \left( \frac{\partial U}{\partial V} \right)_T + P \right] dV$ |
| Adiabatic  | $\delta q = 0$ | $C_V dT = -\left[ \left( \frac{\partial U}{\partial V} \right)_T + P \right] dV$ |

### Define the System

- Ideal gas
  - $U = U(T)$
  - $P = N k_{\text{B}} T / V$
- Monoatomic ideal gas
  - $C_V = \frac{3}{2} N k_{\text{B}}$

### Integrate the First Law

#### Isochoric Process

```{math}
\begin{aligned}
\delta q &= C_V dT \\
\int_{T_1}^{T_2} \delta q &= \int_{T_1}^{T_2} C_V dT \\
q &= \int_{T_1}^{T_2} \frac{3}{2} N k_{\text{B}} dT \\
q &= \frac{3}{2} N k_{\text{B}} \left( T_2 - T_1 \right)
\end{aligned}
```

#### Isothermal Process

```{math}
\begin{aligned}
\delta q &= \left[ \left( \frac{\partial U}{\partial V} \right)_T + P \right] dV \\
\int_{V_1}^{V_2} \delta q &= \int_{V_1}^{V_2} P dV \\
q &= \int_{V_1}^{V_2} \frac{N k_{\text{B}} T}{V} dV \\
q &= N k_{\text{B}} T \ln \left( \frac{V_2}{V_1} \right)
\end{aligned}
```

#### Adiabatic Process

```{math}
\begin{aligned}
C_V dT &= -\left[ \left( \frac{\partial U}{\partial V} \right)_T + P \right] dV \\
\int_{T_1}^{T_2} C_V dT &= -\int_{V_1}^{V_2} P dV \\
\int_{T_1}^{T_2} \frac{C_V}{T} dT &= -\int_{V_1}^{V_2} \frac{N k_{\text{B}}}{V} dV \\
\frac{3}{2} \ln \left( \frac{T_2}{T_1} \right) &= \ln \left( \frac{V_1}{V_2} \right) \\
T_1 V_1^{2/3} &= T_2 V_2^{2/3}
\end{aligned}
```

## Microscopic Interpretation of the First Law

```{math}
\begin{aligned}
U &= \sum_{i = 1}^M p_i E_i \\
dU &= \sum_{i = 1}^M E_i dp_i + \sum_{i = 1}^M p_i dE_i \\
&= \sum_{i = 1}^M E_i dp_i + \sum_{i = 1}^M p_i \left( \frac{\partial E_i}{\partial N} \right)_V dN + \sum_{i = 1}^M p_i \left( \frac{\partial E_i}{\partial V} \right)_N dV \\
&= \sum_{i = 1}^M E_i dp_i + \biggl\langle \left( \frac{\partial E}{\partial V} \right)_N \biggr\rangle dV \\
&= \sum_{i = 1}^M E_i dp_i - P dV
\end{aligned}
```

```{math}
\begin{aligned}
\delta q &= \sum_{i = 1}^M E_i dp_i \\
\delta w &= \sum_{i = 1}^M p_i dE_i
\end{aligned}
```
