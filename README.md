# Global Quasi Linearization (GQL)

The GQL method [1] is mainly introduced and developed for the reduction of chemical kinetics. The generation of GQL reduced chemistry is based on the homogeneous reacting system which can be expressed as ODE system.

The starting point is the ODE system with thermo-kinetic state vector $\psi$:

```math
\frac{\text{d}\psi}{\text{d}t}=F(\psi),
```
where $F(\psi)$ is the vector of source term of the thermo-kinetic state.

The simulation based on the GQL reduced chemistry is to solve the DAE system [2]:

```math
\textbf{M}_s \frac{\text{d}\psi}{\text{d}t}=F(\psi),
```
where $\textbf{M}_s$ is the mass matrix defined as:

```math
\textbf{M}_s = \begin{pmatrix}
                  Z_s ~ Z_f
                  \end{pmatrix} \cdot \begin{pmatrix}
                  \tilde{Z}_s \\
                  \textbf{0}
                  \end{pmatrix},
```

The GQL_RedChem package requires the Matlab coupled with Cantera. 

References:
[1] Bykov, V., V. Gol'Dshtein, and U. Maas. "Simple global reduction technique based on decomposition approach." Combustion Theory and Modelling 12.2 (2008): 389-405.
[2] Yu, U., V. Bykov and U. Maas, "Global quasi-linearization (gql) versus QSSA for a hydrogen–air auto-ignition problem." Physical Chemistry Chemical Physics 20.16 (2018): 10770–10779.
