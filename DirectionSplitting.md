# Solution Algorithm: Fractional-Step Method with Direction Splitting

This document describes the solution algorithm for the **Navier–Stokes–Brinkman equations** using a **Fractional-Step (projection) method** with a *pressure-correction* formulation.  
A key feature of the method is the use of **direction splitting**, applied both to the viscous momentum step and to the pressure-correction step.

The algorithm advances the solution from time level $n$ to $n+1$ through four main stages.

---

## Preliminary Definitions

The following constants arise from the Crank–Nicolson discretization:

$$
\beta = 1 + \frac{\Delta t \nu}{2k}
$$

$$
\gamma = \frac{\Delta t \nu}{2\beta}
$$

---

## Phase 1: Prediction and RHS Assembly

### 1. Pressure Predictor

$$
\rho_*^{n+1/2} = p^{n-1/2} + \phi^{n-1/2}
$$

### 2. RHS Assembly

$$
\hat{g}^{n+1/2}
= f^{n+1/2}
- \nabla \rho_*^{n+1/2}
- \frac{\nu}{2k} u^n
+ \frac{\nu}{2} (\partial_{xx}\eta^n + \partial_{yy}\zeta^n + \partial_{zz}u^n)
$$

---

## Phase 2: Momentum Equation (Velocity Step)

The viscous operator is split as:

$$
(1 - \gamma \nabla^2)
= (1 - \gamma \partial_{xx})(1 - \gamma \partial_{yy})(1 - \gamma \partial_{zz})
$$

### Intermediate Field

$$
\xi^{n+1} = u^n + \frac{\Delta t}{\beta} \hat{g}^{n+1/2}
$$

### X-Direction Solve

$$
(I - \gamma \partial_{xx})(\eta^{n+1} - \eta^n) = \xi^{n+1} - \eta^n
$$

### Y-Direction Solve

$$
(I - \gamma \partial_{yy})(\zeta^{n+1} - \zeta^n) = \eta^{n+1} - \zeta^n
$$

### Z-Direction Solve

$$
(I - \gamma \partial_{zz})(u^{n+1} - u^n) = \zeta^{n+1} - u^n
$$

---

## Phase 3: Pressure Correction Step

Incompressibility:

$$
\nabla \cdot u^{n+1} = 0
$$

Factored operator:

$$
A = -(1 - \partial_{xx})(1 - \partial_{yy})(1 - \partial_{zz})
$$

### X-Direction Solve

$$
(I - \partial_{xx})\psi
= -\frac{1}{\Delta t}\nabla \cdot u^{n+1}
$$

### Y-Direction Solve

$$
(I - \partial_{yy})\varphi = \psi
$$

### Z-Direction Solve

$$
(I - \partial_{zz})\phi^{n+1/2} = \varphi
$$

---

## Phase 4: Pressure Update

$$
p^{n+1/2} = p^{n-1/2} + \phi^{n+1/2}
$$

---

## Conclusion

After this step, the fields:

- $u^{n+1}$
- $p^{n+1/2}$
- $\eta^{n+1}, \zeta^{n+1}, \phi^{n+1/2}$

are available for reuse in the next timestep.
