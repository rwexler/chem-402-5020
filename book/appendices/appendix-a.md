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

# Appendix A. Differential Calculus

## Exact Differentials vs. State Variables

There's a subtle distinction between an exact differential and a state variable that's worth considering. Formally, we say that "a state variable has an exact differential." Alternatively, "the differential of a state variable is exact." For instance, if we compress a balloon infinitesimally ($dV_m$) while slowly moving it toward a heat source ($dT$), the pressure $P$ would change by:

```{math}
dP = \left( \frac{\partial P}{\partial V_m} \right)_T dV_m + \left( \frac{\partial P}{\partial T} \right)_{V_m} dT
```

This illustrates how the value of a state variable $P$ shifts in response to tiny changes in its independent variables. In short, while it's correct to say "$dP$ is exact (not $P$) because $P$ is a state variable," it's more logically precise to note that "$P$ is a state variable because $dP$ is exact."

## Exact vs. Inexact Differentials

In thermodynamics, we encounter both exact and inexact differentials:

* **Exact differentials**: These have equal mixed partial derivatives and integrate to a unique function (for example, $dU$ for the internal energy $U$).

* **Inexact differentials**: These have unequal mixed partial derivatives and do not integrate to a single function (examples include the differentials for heat and work).

The idea of seeking a method to convert an inexact differential into an exact one is a key step in the derivation of entropy. In fact, by applying an integrating factor—in this case, $1/T$—we can transform the inexact differential of heat $\delta q$ into the exact differential $dS = \delta q/T$. This approach was central to the original formulation of entropy.
