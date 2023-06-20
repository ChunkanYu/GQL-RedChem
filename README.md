# Global Quasi Linearization (GQL)

The GQL method is mainly introduced and developed for the reduction of chemical kinetics. The generation of GQL reduced chemistry is based on the homogeneous reacting system which can be expressed as ODE system.

The starting point is the ODE system with thermo-kinetic state vector $\psi$:

```math
\frac{\text{d}\psi}{\text{d}t}=F(\psi),
```
where $F(\psi)$ is the vector of source term of the thermo-kinetic state.

The simulation based on the GQL reduced chemistry is to solve the DAE system:

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
