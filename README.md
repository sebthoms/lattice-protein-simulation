# Lattice Protein Simulation

This project implements a 2D lattice-based protein folding model
using Monte Carlo (Metropolis) dynamics.

Inspired by Computational Physics -	Nicholas J. Giordano, Hisao Nakanish - Chapter 12.1 "Protein Folding".

## Features
- Configurable protein length or sequence
- Lattice constraints
- Energy calculation via interaction matrix
- Monte Carlo sampling

## Requirements
- Python 3.11+
- numpy
- matplotlib
- ffmpeg

## Status
This is an educational / exploratory project.

## Run it

Run `run.py`. Feel free to change the protein configuration in `config.toml` and customize `run.py` to your needs.

## Examples

![[examples/fixedplot_energy_length30_iterations50000_temp2.png]]
![[examples/fixedplot_length_length30_iterations50000_temp2.png]]
![[examples/varplot_energy_length30_iterations1000000_start10_end1.png]]
![[examples/varplot_length_length30_iterations1000000_start10_end1.png]]
![[examples/varplotavg_energy_length30_iterations100000_start10_end1.png]]
![[examples/varplotavg_length_length30_iterations100000_start10_end1.png]]
![Animation example](examples/length100_iterations100000000_temperature2.mp4)
