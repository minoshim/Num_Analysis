[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/minoshim/Num_Analysis/blob/main/Vlasov/Self_g/Py/g_insta/g_insta.ipynb)


# 1D Self-Gravitating Vlasov Simulation

This notebook solves the **1D self-gravitating Vlasov–Poisson system** using a semi-Lagrangian scheme in phase space \((x, v)\).

---

## 📌 Overview

We solve the Vlasov equation:

\[
\frac{\partial f}{\partial t}
+ v \frac{\partial f}{\partial x}
- g(x,t) \frac{\partial f}{\partial v} = 0
\]

with the self-consistent gravitational field:

\[
\frac{\partial g}{\partial x} = -4\pi (\rho - \bar{\rho}), \quad
\rho(x,t) = \int f(x,v,t)\, dv
\]

---

## ⚙️ Numerical Method

- **Operator splitting (Strang splitting)**:
  - Half step: advection in \(x\)
  - Full step: advection in \(v\)
  - Half step: advection in \(x\)

- **Advection scheme**:
  - 3rd-order conservative semi-Lagrangian (CSL3rd)
  - CSL3rd + MUSCL-type reconstruction with limiter (CSLMSL)

- **Boundary conditions**:
  - Periodic in \(x\)
  - Free (no update) in \(v\)

---

## 🧮 Simulation Setup

### Grid

| Parameter | Meaning |
|----------|--------|
| `xmesh` | Number of spatial cells |
| `vmesh` | Number of velocity cells |
| `xoff`, `voff` | Ghost cell width |

Total grid:
  - nx = xmesh + 2xoff
  - nv = vmesh + 2voff


---

### Domain

- Spatial domain: \(x \in [-L/2, L/2]\)
- Velocity domain: \(v \in [-V, V]\)

---

### Initial Condition

Maxwellian in velocity with sinusoidal perturbation in space:

\[
f(x,v,0) = f_v(v)\left[1 + A \cos(kx)\right]
\]

- `amp`: perturbation amplitude
- `kk`: perturbation wavenumber
- `vs`: thermal velocity

---

### Time Step

\[
dt = \text{CFL} \cdot \frac{dx}{v_{\max}}
\]

---

## 🚀 How to Run

1. Open the notebook in Colab
2. Modify parameters at the top:
3. Run all cells
