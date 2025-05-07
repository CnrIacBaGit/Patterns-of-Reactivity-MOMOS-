# Patterns-of-Reactivity-MOMOS-
Generate patterns in the asymptotically stable and reactive region of the parameter space.

# MOMOS Model Simulation

## Overview

This MATLAB script simulates the MOMOS reaction-diffusion system with chemotaxis. The system consists of two coupled PDEs for variables u and v, involving diffusion, chemotactic interaction, and nonlinear reaction kinetics.

The user is prompted at runtime to choose one of four predefined parameter sets, each corresponding to a different regime in the parameter space. The corresponding initial conditions are automatically applied.

---

## Model Equations

The system being solved is:

    u_t = Du * ∆u - div(beta * h(u) * grad(v)) + f(u, v)
    v_t = Dv * ∆v + g(u, v)

with:
- h(u) = u
- f(u,v) = -k1 * u - q * u * |u| + k2 * v
- g(u,v) = k1 * u - k2 * v + c

---

## How to Run

1. Open the script in MATLAB.
2. Run the file.
3. When prompted, enter a number between 1 and 4 to select a case.
4. The simulation will execute and display the final results at T = 500.

---

## Parameter Sets

| Case | q         | beta       | sigma | Initial Conditions           |
|------|-----------|------------|--------|------------------------------|
| 1    | 0.0433    | 0.806      | 10     | Random perturbation          |
| 2    | 0.061122  | 1.01668    | 12     | Top band with constant value |
| 3    | 0.0196639 | 0.474095   | 12     | Random perturbation          |
| 4    | 0.0804361 | 1.23535    | 7      | Top band with constant value |

- Random perturbation: small noise around equilibrium.
- Top band: top 30% of the domain set to a constant high value.

---

## Output

Two figures are shown at the end:
- `u at T = 500`: distribution of u
- `v at T = 500`: distribution of v

Both figures use interpolated shading, a jet colormap, and font size 14 with bold labels.

---

## Requirements

- MATLAB (R2018a or later recommended)
- No external toolboxes required

---

## Notes

- The spatial domain is square with periodic boundary conditions.
- Finite difference discretization is used for space; Kronecker products for 2D operators.
- Time integration is performed using a Implicit Symplectic Euler scheme with LU decomposition.

---

## Contact

For any questions or improvements, feel free to contact the author.
