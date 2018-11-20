# Function help

## 1. lagTensegrityDynamics
Main ODE function generating the algebraic differential equations to be solved using ODE45 and the novel constraint correction method.

## 2. ConstraintCorrection
Function to minimise constraint violations occurring due to numerical integration. Called after a successful step from ODE45.m

## 3. tensegConstruct
Function to create the essential parameters that uniquely define the tensegrity structure.

## 4. tensegGenMat
Function to generate matrices under the Lagrangian framework described in the [preprint].

## 5. tensegEq
Function to find the equilibrium force densities in the strings at t=0, given initial position and minimum force densities in the strings.
The initial nodal coordinates (q) are assumed to be given by the user.
In equation (31) from the [preprint], the first term on the left (\ddot{q}) is taken to be 0.
The resulting equation is solved for \sigma and \lambda while constraining the force densities (\sigma) to be positive and above the lower bound set by the user.
The equilibrium force densities are used to find the new rest lengths of the springs from equation 25, and subsequently, cable energies (equation 26) to be used in the Lagrangian dynamics.


## 6. tensegSimTime
Function to generate simulation time and output interval time-step for simulation and animation.

## 7. tensegSim
Function to simulate the tensegrity structure's dynamics in the prescribed environment.

## 8. plot_configuration
Function to plot initial configuration of the tensegrity structure.

## 9. plotMotion
Function to plot trajectories of all nodes in the tensegrity structure after solution has been obtained.

## 10. plotConstr
Function to plot constraint violations over time.

## 11. animateTenseg
Function to animate the simulated dynamics of the tensegrity structure.

[preprint]: https://www.researchgate.net/publication/328676032_A_Lagrangian_Formulation_for_Constrained_Multibody_Dynamics_in_Tensegrity_Systems
