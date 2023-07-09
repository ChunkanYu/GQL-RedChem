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

The "Davis_and_Skodje_Model" includes one simple case example (coming soon).

In "GQL_RedChem" package, one can find the main code for the generation of GQL reduced chemistry. They require the Matlab coupled with Cantera.

How to use?
1) main_0dSimulation.m: calculate the detailed chemistry, GQL reduced chemistry and QSSA reduced chemistry.
2) main_GQL_RedChem.m: code to find out the GQL reduced chemistry, which will be stored in GQL_Ms.mat. Note that there are many candidates stored in the GQL_Ms.mat. One needs to select and test them in main_0dSimulation for his validation for a wider range of application.

References and List of all scholarly publications enabled by the software till 09th July 2023:

[1] V. Bykov., V. Gol'Dshtein, and U. Maas. "Simple global reduction technique based on decomposition approach", Combustion Theory and Modelling 12.2 (2008): 389-405.
[2] C. Yu, V. Bykov and U. Maas, "Global quasi-linearization (gql) versus QSSA for a hydrogen–air auto-ignition problem", Physical Chemistry Chemical Physics 20.16 (2018): 10770–10779.
[3] V. Bykov, C. Yu, V. Gol’dshtein, U. Maas, "Model reduction and mechanism comparison of hydrogen/oxygen auto-ignition", Proceedings of the Combustion Institute 37.1 (2019): 781 – 787.
[4] C. Yu, V. Bykov, U. Maas, "Coupling of simplified chemistry with mixing processes in pdf simulations of turbulent flames", Proceedings of the Combustion Institute 37.2 (2019): 2183–2190.
[5] C. Yu, F. Minuzzi, V. Bykov, U. Maas, "Methane/air auto-ignition based on global quasi-linearization (GQL) and directed relation graph (DRG): Implementation and comparison", Combustion Science and Technology 192.9 (2020): 1802–1824.
[6] V. Bykov, S. Shashidharan, E. Berszany, V. Gubernov, U. Maas, "Model reduction of rich premixed hydrogen/air oscillatory flames by global quasi-inearization (GQL)", Combustion Science and Technology 194.12 (2022): 2377–2394.
