# Numerical Rocket Simulation

Computational project for the simulation of a rocket launch using physical models
and numerical methods, implemented in Python.

## Description
The project addresses the simulation of a rocket’s motion, taking into account mass variation
during combustion and the action of different external forces.  
The model is constructed progressively, starting from simple approximations and
adding more realistic physical effects.

Both single-stage and multi-stage rockets are studied, analyzing the behavior of the
trajectory, velocity, and fuel consumption.

## Features
- Variable-mass models
- Constant gravity and gravity dependent on height and latitude
- Aerodynamic drag with altitude-dependent air density
- Numerical solution of differential equations using different methods
- Comparison between physical models and integration schemes
- Single-stage and multi-stage rocket simulations

## Repository Structure
- `FinalProject.py`: main simulation and results visualization script
- `GlobalFunctionsFin.py`: implementation of physical models and auxiliary functions
- `Rocket_MNAF.pdf`: detailed documentation of the approach and results

## Requirements
- Python 3.x
- NumPy
- SciPy
- Matplotlib

## Usage
Run the main script "FinalProject.py"

## Authors
Ignacio González Álvarez, Miguel Durán Vera, David Villa Blanco.

This repository is for academic purposes.
