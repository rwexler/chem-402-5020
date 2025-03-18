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

```{glossary}
Quasi-static process
  A process carried out *infinitesimally slowly*, ensuring the system remains *near equilibrium* at each step. Any deviations from equilibrium are so small (infinitesimal) that each intermediate state is well-defined thermodynamically.

Reversible process
  An *idealized* process carried out *quasistatically* and with *no dissipative effects*, so that it can be undone without leaving any net change in the system or surroundings.
```

## How to Apply the First Law of Thermodynamics

```{admonition} 1. Choose Two Independent Variables
:class: dropdown

- **How to Choose**:
  - *Convenience*: What experiments can you set up most easily with the equipment at hand?
  - *Necessity*: Which variables are feasible to control or measure with the facilities available?
  - *Curiosity*: Which variables will yield the most interesting or relevant results for your study?
- Generally, two independent variables are sufficient for many thermodynamic experiments; controlling more can become unwieldy.
```

````{admonition} 2. Rewrite the First Law
:class: dropdown

- Express $\delta q$ (the infinitesimal heat transfer) in terms of the differentials of your chosen independent variables, $X_1$ and $X_2$.

```{math}
\delta q = 
\left(\frac{\partial U}{\partial X_1}\right)_{X_2} dX_1 + 
\left(\frac{\partial U}{\partial X_2}\right)_{X_1} dX_2 
+ \dots
```

- If work is done, include an appropriate term for $\delta w$. For example, if $P$–$V$ work is relevant, $\delta w = -P\, dV$.
````

```{admonition} 3. Apply Constraints
:class: dropdown

- **Holding Variables Constant**:
  - *Isobaric*: $dP = 0$
  - *Isochoric*: $dV = 0$
  - *Isothermal*: $dT = 0$
- **Adiabatic Boundary**:
  - If the boundary is thermally insulating, $\delta q = 0$.
  - Note that “adiabatic” means no heat transfer; this does not necessarily imply the system is completely isolated (it can still do work, for example).
```

```{admonition} 4. Define the System
:class: dropdown

- Specify the system’s equation of state if known (e.g., $PV = nRT$ for an ideal gas).
- The equation of state helps you calculate partial derivatives (e.g., $\left(\frac{\partial U}{\partial V}\right)_T$) that appear in your expression for $\delta q$.
```

```{admonition} 5. Integrate the First Law
:class: dropdown

- By integrating, you eliminate the differentials to get explicit relationships. 
- In real experiments, processes are rarely perfectly quasi-static or reversible, but the differential approach can still help interpret results.
- The final expression links the heat exchanged $(q)$ to the changes in your chosen independent variables along the process path.
```

```{admonition} So What?
:class: tip

By carefully selecting your independent variables, rewriting the First Law to account for them, and integrating under specified constraints, you can predict how energy conservation will impact measurable quantities such as heat absorbed or released. Different constraints (isobaric, isochoric, isothermal, or adiabatic) reveal how changing one variable can affect others—and the corresponding energy flow.
```

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
