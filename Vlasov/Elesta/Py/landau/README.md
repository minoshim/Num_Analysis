[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/minoshim/Num_Analysis/blob/main/Vlasov/Elesta/Py/landau/landau.ipynb)

# 1D electrostatic Vlasov Simulation

This notebook solves the **1D electrostatic Vlasov–Poisson system** using a semi-Lagrangian scheme in phase space \((x, v)\).

---

## 📌 Overview

We solve the Vlasov equation for the phase space distiribution f(x,v):

df/dt + vdf/dx qEdf/dv = 0,

with the self-consistent electric field E(x):

dE/dx = q (rho - rho_0), 

where
  - rho: number density (integration of f in velocity space)
  - rho_0: mean number density

---

## ⚙️ Numerical Method

- **Operator splitting (Strang splitting)**:
  - Half step: advection in \(x\)
  - Full step: advection in \(v\)
  - Half step: advection in \(x\)

- **Advection scheme**:
  - 3rd-order conservative semi-Lagrangian (CSL3rd), or
  - CSL3rd + MUSCL-type flux limiter (CSLMSL)

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

This example simulates the two-stream instability:

f(x,v,0) = 0.5(f_v+(v) + f_v-(v)) [1 + A cos(kx)]

where f_v+(v) and f_v-(v) are the shifted gaussian (positive and negative drift velocity)
- `amp`: perturbation amplitude
- `kk`: perturbation wavenumber
- `vd`: drift velocity
- `vs`: thermal velocity
  
You may simulate other tests such as linear and nonlinear Landau damping using different initial condition.

---


## 🚀 How to Run

1. Open the notebook in Colab: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/minoshim/Num_Analysis/blob/main/Vlasov/Elesta/Py/landau/landau.ipynb)
2. Modify parameters at the top:
3. Run all cells
