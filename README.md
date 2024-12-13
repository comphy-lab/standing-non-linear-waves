# Standing Non-Linear Waves

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This repository contains simulation code for studying the dynamics of standing non-linear waves using the Basilisk flow solver. The code implements both analytical Stokes wave solutions and experimental best-fit initial conditions to investigate wave dynamics with adaptive mesh refinement.

## Overview

The project focuses on simulating standing waves using a two-phase flow solver with surface tension. It explores:
- Non-linear wave dynamics using Stokes wave theory
- Adaptive mesh refinement for accurate interface tracking
- Comparison between analytical and experimental initial conditions

## Installation and Setup

To ensure you have the necessary tools and a fresh Basilisk installation, use the provided script:

```bash
./reset_install_requirements.sh
```

### Function
This script checks for Basilisk installation and compiles it if not present.

### OS Compatibility
Designed for macOS and Linux systems.

### Dependencies
- Basilisk C (automatically installed)
- C compiler (gcc/clang)
- Make build system

### Environment Setup
After running the script, a `.project_config` file is created, setting `BASILISK` and `PATH` automatically.

For a fresh installation:
```bash
./reset_install_requirements.sh --hard
```

## Running the Code

### Using Makefile (Recommended)

1. Navigate to the test cases directory:
```bash
cd testCases
```

2. Compile and run:
```bash
make StokesStandingWaves
```

### Direct Compilation

```bash
qcc -O2 -Wall -disable-dimensions StokesStandingWaves.c -o StokesStandingWaves -lm
```

### Execution

The program accepts the following parameters:
```bash
./StokesStandingWaves maxLevel We Ohd tf Ohs De Ec tmax
```

where:
- `maxLevel`: Maximum refinement level for adaptive mesh
- `We`: Weber number
- `Ohd`: Ohnesorge number for the drop
- `tf`: Shell thickness
- `Ohs`: Shell Ohnesorge number
- `De`: Deborah number
- `Ec`: Elasto-capillary number
- `tmax`: Maximum simulation time

## Features

- Two-phase flow simulation with surface tension
- Adaptive mesh refinement for interface tracking
- Support for both analytical and experimental initial conditions
- Automatic error control and mesh adaptation
- Comprehensive logging and data output

## Technical Details

### Key Components
- Volume of Fluid (VoF) method for interface tracking
- Height function method for curvature calculation
- Adaptive mesh refinement with error control
- Stokes wave analytical solutions

### Simulation Parameters
- Domain size: 2.0 × 2.0
- Error tolerances:
  - VOF: 1e-3
  - Curvature: 1e-6
  - Velocity: 1e-3

## Project Team

- Vatsal Sanjay (University of Twente)
  - Email: v.sanjay@utwente.nl

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this code in your research, please cite:

```bibtex
@software{sanjay2024standing,
    author = {Vatsal Sanjay},
    title = {Standing Non-Linear Waves},
    year = {2024},
    url = {https://github.com/VatsalSy/standing-non-linear-waves}
}