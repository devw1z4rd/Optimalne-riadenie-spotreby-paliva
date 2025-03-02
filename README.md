# Optimal Fuel Consumption Control

## Abstract
This MATLAB script, `solve_fuel_minimization_corrected.m`, addresses the problem of optimal fuel consumption control over a finite time horizon (`N = 30`). The approach leverages linear programming (LP) to minimize an auxiliary decision variable, `z_t`, which encapsulates a piecewise-defined cost function `f(u(t)) = max(|u|, 2|u| - 1)`. The formulation ensures compliance with the system's dynamic constraints while achieving a specified terminal state.

## Problem Statement
- **System Dynamics:** The evolution of the state vector is governed by the linear discrete-time equation with transition matrix `A` and control influence vector `b`.
- **Terminal Constraint:** The state at the final time step must satisfy `x(N) = [7; 2; -6]`.
- **Optimization Objective:** The objective function is constructed to minimize the cumulative sum of `z_t` over the control horizon, ensuring adherence to predefined constraints that regulate control effort and state evolution.

## Computational Approach
1. Reformulates the fuel consumption minimization problem as a linear program with non-negative auxiliary variables `u_t^+`, `u_t^-`, and `z_t`.
2. Constructs equality constraints to enforce system dynamics and final state conditions.
3. Imposes inequality constraints to enforce the structure of the cost function.
4. Utilizes MATLAB’s `linprog()` function under the `dual-simplex` algorithm for efficient solution computation.
5. Post-processes the results to extract optimal control inputs `u(t)`, auxiliary variables `z(t)`, and the corresponding state trajectory `x(t)`.
6. Provides a visualization of the computed control sequence and state evolution over time.

## Execution Instructions
To execute the script in MATLAB, run:
```matlab
solve_fuel_minimization_corrected
```

## Prerequisites
- MATLAB with the `Optimization Toolbox` installed.

## Acknowledgments
This implementation was developed as part of a study on fuel-efficient control strategies within constrained dynamic systems.
