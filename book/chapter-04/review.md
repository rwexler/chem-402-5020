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

# Review

## 1. Checklist of Key Concepts

1. **Entropy as a State Function**
   - **Definition**: For a reversible process, the change in entropy is

   ```{math}
   dS \;=\;\frac{\delta q_\text{rev}}{T}
   ```

   - **Key Property**: $S$ depends only on the initial and final states—**not** on the path taken.
   - **Physical Meaning**: Measures dispersal of energy; quantifies how heat “spreads out” in a system.

2. **Spontaneity and the Second Law**
   - **Clausius Statement**: Heat does not spontaneously flow from cold to hot.
   - **Entropy Criterion**:
     - For an **isolated** system, $\Delta S_\text{tot} \ge 0$.
     - Equality holds only for reversible processes; strict inequality for irreversible ones.
   - **Macroscopic Implication**: Determines the “arrow of time”—which processes can occur spontaneously.

3. **Carnot Cycle & Maximum Efficiency**
   - **Reversible Engine**: Comprised of two isotherms $\left( T_h, T_c \right)$ and two adiabats.
   - **Efficiency**:

   ```{math}
   \eta_\text{Carnot} \;=\; 1 - \frac{T_c}{T_h}
   \quad\Longrightarrow\quad
   \text{This is the upper limit for any heat engine.}
   ```

   - **Significance**: Links thermodynamic reversibility to practical limits on converting heat to work.

4. **Microscopic Origin of Entropy**
   - **Boltzmann’s Formula**:

   ```{math}
   S \;=\; k_B \ln \Omega
   ```

   where $\Omega$ is the number of accessible microstates at given $U, V, N$.
   - **Statistical Interpretation**: Entropy quantifies uncertainty about which microstate the system occupies.
   - **Approach to Equilibrium**: Systems evolve towards macrostates with higher $\Omega$; thus $S$ increases until equilibrium.

---

## 2. Checklist of Most Important Equations

| **Equation** | **Variable Definitions** | **Meaning & Use Cases** | **Caveats/Approximations** |
|:-|:-|:-|:-|
| $dS = \frac{\delta q_\text{rev}}{T}$ | $\delta q_\text{rev}$: infinitesimal heat in a reversible process<br>$T$: absolute temperature | Fundamental relation defining entropy. Use to compute $\Delta S$ for isothermal, isobaric, mixing, etc. processes. | Valid **only** for reversible paths; for irreversible processes, $\Delta S > \int \delta q/T$. |
| $dU = T\,dS - P\,dV$ | $U$: internal energy<br>$P$: pressure<br>$V$: volume | **First law** written in natural variables $(S, V)$. Use to derive Maxwell relations and fundamental thermodynamic identities. | Assumes $U=U(S,V)$ is a well-defined state function; neglects non-PV work. |
| $\eta_\text{Carnot} = 1 - \dfrac{T_c}{T_h}$ | $T_h$, $T_c$: temperatures of hot and cold reservoirs | Sets the **maximum** efficiency for any heat engine operating between $T_h$ and $T_c$. Use in engine design and in defining thermodynamic reversibility. | Requires fully reversible (Carnot) cycle; real engines incur friction, heat leaks, and irreversibilities. |
| $\Delta S_\text{tot} \ge 0$ | $\Delta S_\text{tot}$: entropy change of system + surroundings | General statement of the **Second Law** for isolated systems. Predicts direction of spontaneous processes. | Equality only for reversible processes; in open or non-isolated systems, must include entropy change of environment. |
| $S = k_B \ln \Omega$ | $k_B$: Boltzmann constant<br>$\Omega$: number of microstates | Bridges macroscopic entropy with molecular statistics. Use to calculate $S$ from combinatorial counts or partition functions. | Applies to **isolated** systems with well-defined $\Omega$. Counting microstates can be intractable for large systems. |
| $p_i = \dfrac{e^{-\beta E_i}}{Q}$; $\beta = 1/(k_B T)$ | $E_i$: energy of microstate $i$<br>$Q$: canonical partition function | Probability of occupying microstate $i$ in the canonical ensemble. Use in statistical mechanics to compute ensemble averages. | Assumes equilibrium with a heat bath at fixed $T$; neglects quantum degeneracy unless included in $E_i$. |
| $S = \frac{U}{T} + k_B \ln Q$ | $Q = \sum_i e^{-\beta E_i}$: partition function | Expresses entropy in terms of the partition function and internal energy. Use to derive free energies, heat capacities, and response functions. | Derived under canonical ensemble assumptions. Valid when $U$ and $Q$ are known or can be approximated. |
