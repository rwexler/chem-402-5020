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

# 3.1. Conservation of Energy

> *“The total energy of the world, kinetic plus potential, is a constant when we look closely enough… [If we see energy not conserved, then] this is due to a lack of appreciation of what it is that we see.”*  
> — Richard Feynman

## Overview

In this section, we introduce the principle of **conservation of energy** and the **First Law of Thermodynamics**, which unifies the ideas of heat and work into one fundamental statement.

## Conservation of Mechanical Energy

Between 1732 and 1736, Bernoulli and Euler combined the discoveries of Newton (laws of motion) and Leibniz (the connection between weight × vertical displacement and weight × velocity squared) into an early form of the law of conservation of mechanical energy. A simple example is the interchange between potential and kinetic energy:

```{math}
\frac{1}{2} m v^2 = m g h.
```

When an object (mass $m$) of height $h$ is dropped, its gravitational potential energy $mgh$ is converted into kinetic energy $\tfrac{1}{2}mv^2$, assuming negligible air resistance.

````{admonition} Apollo 15
:class: tip

During the Apollo 15 mission, astronaut David Scott famously dropped a hammer and a feather on the Moon. With negligible air resistance, both hit the ground simultaneously. This illustrates that in free fall, all objects (no matter their mass) accelerate at the same rate, and thus gain the same velocity over the same drop distance:

```{math}
\frac{1}{2}\,\cancel{m} v^2 = \cancel{m} g h
\quad\Rightarrow\quad
v = \sqrt{2 g h}.
```

```{image} https://upload.wikimedia.org/wikipedia/commons/6/62/AS15-88-11866_-_Apollo_15_flag%2C_rover%2C_LM%2C_Irwin_-_restoration1.jpg
:width: 50%
:align: center
```
````

## Mechanical Equivalent of Heat

```{figure} https://upload.wikimedia.org/wikipedia/commons/c/c3/Joule%27s_Apparatus_%28Harper%27s_Scan%29.png
---
name: joule-apparatus
width: 50%
align: center
---
Joule's apparatus for measuring the mechanical equivalent of heat.
```

In 1847, James Prescott Joule measured how mechanical work converts into heat. He famously found that dropping a total of $778$ lb$\cdot$ft of weight (e.g., by turning a paddle in water) raised the temperature of $1$ lb of water by $1\,^\circ\text{F}$. Converting to SI units:

```{math}
\text{MEH} 
= \frac{778 \;\text{lb}\cdot\text{ft}}{1\;\text{lb}\,^\circ\text{F}}
\approx \frac{1{,}055 \;\text{J}}{(453.6 \;\text{g water}) \cdot (5/9\,^\circ\text{C})}
\approx 4.18 \;\text{J}\,\text{g}^{-1}\,{^\circ\text{C}}^{-1}.
```

Here, $4.18\,\text{J}\,\text{g}^{-1}\,{^\circ\text{C}}^{-1}$ is the **specific heat** of water—i.e., the heat capacity per gram.

## First Law of Thermodynamics

Joule’s work revealed that heat and mechanical energy are interconvertible. To include heat ($q$) and work ($w$) in one statement of energy conservation, we use:

```{math}
\Delta U = q + w,
```

where

- $\Delta U$ is the **change in internal energy**,
- $q$ is the heat **absorbed by** the system,
- $w$ is the work **done on** the system.

### Differential Form

In differential form,

```{math}
dU = \delta q + \delta w 
\quad\Longleftrightarrow\quad
\boxed{\delta q = dU - \delta w.}
```

The notation $\delta$ is used to remind us that $q$ and $w$ are **path functions**, not state functions.

## Types of Work

In Module 1.1, we defined work as “energy transferred when a force acts over a distance.” Mathematically:

```{math}
\delta w = \vec{F}\cdot d\vec{r} 
= (F_x,\,F_y,\,F_z)\cdot(dx,\,dy,\,dz) 
= F_x\,dx + F_y\,dy + F_z\,dz.
```

Many physical processes fit a “generalized force” $\times$ “generalized displacement” pattern:

| Generalized “Force”   | Generalized “Displacement” | $\delta w$            | Example                         |
|:---------------------:|:--------------------------:|:-----------------------:|:--------------------------------:|
| Mechanical $F$       | $x$                     | $F \,dx$              | Lifting a weight                |
| Linear tension $k$   | $l=x - x_0$             | $k \,dl$              | Stretching a spring             |
| Surface tension $\gamma$ | $A$                 | $\pm\,\gamma \, dA$   | Blowing a soap bubble           |
| Pressure $P$         | $V$                     | $-P\,dV$              | Compressing a gas               |
| Chemical $\mu$       | $N$ or $n$            | $\mu\,dN$             | Forming a molecule              |
| Electrical $E$       | Charge $q_{\text{el}}$  | $E\,dq_{\text{el}}$   | Moving an electric charge       |

**Sign convention** in chemistry:  

- $w>0$ when work is **done on** the system (system’s energy goes up).  
- $w<0$ when work is **done by** the system (system’s energy goes down).

---

### Example 1. Expanding Against Constant Pressure

Calculate the work necessary to expand an ideal gas from $20\,\mathrm{L}$ to $85\,\mathrm{L}$ against a constant external pressure of $2.5\,\mathrm{bar}$.  

````{admonition} Solution
:class: dropdown

1. **Convert Units to SI**  
   ```{math}
   P_{\mathrm{ext}} = 2.5\,\mathrm{bar} = 2.5\times10^5\,\mathrm{Pa},\quad
   V_i = 20\,\mathrm{L} = 20\times10^{-3}\,\mathrm{m^3},\quad
   V_f = 85\,\mathrm{L} = 85\times10^{-3}\,\mathrm{m^3}.
   ```

2. **Set Up the Work Integral**  
   For a constant external pressure,
   ```{math}
   \delta w = -P_{\mathrm{ext}}\,dV 
   \quad \Longrightarrow\quad
   w = -\int_{V_i}^{V_f} P_{\mathrm{ext}}\,dV.
   ```

3. **Perform the Integration**  
   ```{math}
   w 
   = -P_{\mathrm{ext}}\int_{V_i}^{V_f} dV
   = -P_{\mathrm{ext}}\bigl(V_f - V_i\bigr).
   ```

4. **Numerical Result**  
   ```{math}
   w 
   = -\bigl(2.5\times10^5\,\mathrm{Pa}\bigr)
     \Bigl[\bigl(85\times10^{-3}\bigr) - \bigl(20\times10^{-3}\bigr)\Bigr]\mathrm{m^3}
   = -16{,}250\,\mathrm{J}
   = \boxed{-16\,\mathrm{kJ}}.
   ```

5. **Interpretation**  
   - The **negative sign** indicates the system (gas) does $16\,\mathrm{kJ}$ of work **on the surroundings**.
   - According to the convention $\Delta U = q + w$, if no heat is exchanged ($q=0$), the gas would lose $16\,\mathrm{kJ}$ in internal energy because it expands.
````

---

### Example 2. Expanding a Water Bubble

Calculate the work necessary to expand a water “bubble” (two surfaces) with surface tension $\gamma = 72\,\mathrm{J/m^2}$ from a radius of $1\,\mathrm{cm}$ to $3.25\,\mathrm{cm}$.

````{admonition} Solution
:class: dropdown

1. **Convert Units**  
   ```{math}
   \gamma = 72\,\mathrm{J/m^2}, \quad
   r_i = 1\times10^{-2}\,\mathrm{m}, \quad
   r_f = 3.25\times10^{-2}\,\mathrm{m}.
   ```

2. **Identify the Surface Area Change**  
   - For a **single** spherical surface, $A = 4\pi r^2$, so $dA = 8\pi r\,dr$.  
   - A bubble has **two** surfaces (inner and outer), so its total area is $2 \times 4\pi r^2 = 8\pi r^2$. Hence,
     ```{math}
     \text{(bubble area)} = 8\pi r^2 
     \quad\Longrightarrow\quad
     dA = 16\pi r\,dr.
     ```

3. **Write the Work Expression**  
   Since $\delta w = -\gamma\,dA$ (negative sign for expansion work done **by** the system):
   ```{math}
   w = - \int_{A_i}^{A_f} \gamma \,dA 
     = - \gamma \bigl(A_f - A_i\bigr).
   ```

4. **Perform the Integral**  
   Substituting $A = 8\pi r^2$ for the bubble:
   ```{math}
   \Delta A = A_f - A_i
   = 8\pi \bigl(r_f^2 - r_i^2\bigr).
   ```
   Thus,
   ```{math}
   w = -\,\gamma\,\Delta A 
      = -\,\gamma \,8\pi\bigl(r_f^2 - r_i^2\bigr).
   ```

5. **Calculate Numerically**  
   ```{math}
   w 
   = -\,72\,\mathrm{J/m^2}
     \,\times\,8\pi
     \Bigl[\bigl(3.25\times10^{-2}\bigr)^2 
           - \bigl(1\times10^{-2}\bigr)^2\Bigr]\mathrm{m^2}.
   ```
   Evaluate numerically (you should find):
   ```{math}
   w \approx \boxed{-1.7\,\mathrm{J}}.
   ```

6. **Interpretation**  
   - Expanding the bubble (from smaller to larger radius) **does work on** the surroundings, giving $w<0$ from the system’s perspective.
   - If no heat is supplied, the system’s internal energy would decrease by $1.7\,\mathrm{J}$.
````

---

### Example 3. Stretching a Hookean Fiber

Calculate the work required to stretch a fiber obeying Hooke’s law, with $k=100\,\mathrm{N/cm}$, by $0.15\,\mathrm{cm}$.

````{admonition} Solution
:class: dropdown

1. **Convert Units**  
   ```{math}
   k = 100\,\frac{\mathrm{N}}{\mathrm{cm}}
     = 100\times10^2\,\frac{\mathrm{N}}{\mathrm{m}}
     = 10^4\,\frac{\mathrm{N}}{\mathrm{m}}, 
     \quad
   l = 0.15\,\mathrm{cm} = 0.15\times10^{-2}\,\mathrm{m}.
   ```

2. **Write the Work Expression**  
   Hooke’s law for tension: $F = k\,\Delta l$. The infinitesimal work is:
   ```{math}
   \delta w = F\,dx = k(x - x_0)\,dx.
   ```

3. **Integrate**  
   If we stretch from $x_0$ to $x = x_0 + l$:
   ```{math}
   w 
   = \int_{x_0}^{x_0 + l} k\,(x - x_0)\,dx
   = \frac{1}{2} \,k\,(l)^2.
   ```

4. **Substitute Numbers**  
   ```{math}
   w 
   = \tfrac{1}{2}\,(10^4)\,
     \Bigl(0.15\times10^{-2}\,\mathrm{m}\Bigr)^2 
   = 0.5\times10^4
     \times\bigl(0.15\times10^{-2}\bigr)^2\,\mathrm{J}.
   ```
   Numerically,
   ```{math}
   w \approx 1.1\times10^{-2}\,\mathrm{J}.
   ```

5. **Interpretation**  
   - This is the energy required (work done **on** the system) to stretch the fiber by $0.15\,\mathrm{cm}$.
   - The positive sign indicates the system absorbs energy (an external force pulls on the fiber).
````
