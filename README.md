# TrebuchetDynamicsSim

## Overview
TrebuchetDynamicsSim is a physics-based simulation of a trebuchet mechanism, developed as part of ME EN 534 - **Dynamics of Mechanical Systems** at Brigham Young University. The project models the coupled motion of the trebuchet components using classical dynamics and numerically simulates the system response over time.

The simulator demonstrates how mass distribution, geometry, and initial conditions influence the motion and projectile launch behavior.

---

## Features
- Dynamic modeling of a multi-body mechanical system
- Time-domain numerical simulation of trebuchet motion
- Visualization of system kinematics
- Parameterized model for exploring design changes

---

## Modeling Approach
The system is modeled using Lagranges Method.  The equations of motion are derived based on:
- Mass and inertia properties of the trebuchet components
- Gravitational forces and constraints
- Kinematic relationships between rotating links

The resulting equations are integrated numerically to simulate the time evolution of the system.

---

## Results
The simulation successfully reproduces realistic trebuchet motion and illustrates the sensitivity of projectile launch behavior to system parameters such as counterweight mass and arm geometry.

---

## How to Run
1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/TrebuchetDynamicsSim.git
   ```
2. Run Trebuchet_Simulation.m
