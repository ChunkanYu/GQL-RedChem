# Davis and Skodje Model

This project consists of the MatLAB code to solve the Davis-Skodje Model proposed in [1].

```math
\frac{\text{d}\psi_1}{\text{d}t}=-\psi_1
```

```math
\frac{\text{d}\psi_2}{\text{d}t}=- \gamma \psi_2 + \frac{(\gamma-1)\psi_1 + \gamma \psi_1^2}{(1 + \psi_1)^2}
```

where $\gamma>1$ gives the stiffness of the ODE system. 

The following methods for the reduction of the chemical kinetics are used:
* the Quasi-steady state approximation (QSSA) [2]
* the Intrinsic Low-Dimensional Manifold (ILDM) [3]
* the Global Quasi-Linearization (GQL) [4]

The main.m file provides:
* the detailed solution
* the QSSA solution
* the ILDM solution
* the GQL solution with 0th order approximation



References:

[1] Davis, Michael J., and Rex T. Skodje. "Geometric investigation of low-dimensional manifolds in systems approaching equilibrium." The Journal of chemical physics 111.3 (1999): 859-874.

[2] Bodenstein, Max. "Eine theorie der photochemischen reaktionsgeschwindigkeiten." Zeitschrift für physikalische Chemie 85.1 (1913): 329-397.

[3] Maas, Ulrich, and Stephen B. Pope. "Simplifying chemical kinetics: intrinsic low-dimensional manifolds in composition space." Combustion and flame 88.3-4 (1992): 239-264.

[4] Bykov, V., V. Gol'Dshtein, and U. Maas. "Simple global reduction technique based on decomposition approach." Combustion Theory and Modelling 12.2 (2008): 389-405.

