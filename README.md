# Standing Non-Linear Waves

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub issues](https://img.shields.io/github/issues/comphy-lab/standing-non-linear-waves)](https://github.com/comphy-lab/standing-non-linear-waves/issues)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14551990.svg)](https://doi.org/10.5281/zenodo.14551990)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/comphy-lab/standing-non-linear-waves)](https://github.com/comphy-lab/standing-non-linear-waves/releases)
[![GitHub last commit](https://img.shields.io/github/last-commit/comphy-lab/standing-non-linear-waves)](https://github.com/comphy-lab/standing-non-linear-waves/commits)
[![Language](https://img.shields.io/github/languages/top/comphy-lab/standing-non-linear-waves)](https://github.com/comphy-lab/standing-non-linear-waves)

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
./StokesStandingWaves maxLevel Ga Bo A0 ORDER tmax
```

where:
- `maxLevel`: Maximum refinement level for adaptive mesh (default: 7)
- `Ga`: Gallileo number (ratio of gravitational to viscous forces)
- `Bo`: Bond number (ratio of gravitational to surface tension forces)
- `A0`: Amplitude of the standing wave
- `ORDER`: Order of the initial condition (0-8, or -1 for best fit)
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

### Physical Parameters
- Density ratio (rho1/rho2): 1000 (water-air like)
- Viscosity ratio (mu2/mu1): 0.01 (water-air like)
- Surface tension coefficient: 1.0/Bo
- Gravity: -1.0 (dimensionless)

### Simulation Parameters
- Domain size: 2.0 × 2.0
- Error tolerances:
  - VOF: 1e-3 (interface tracking)
  - Curvature: 1e-6 (height function method)
  - Velocity: 1e-3 (adjust based on Oh number)
- Grid resolution: 128 points per unit length (in post-processing)

## Post-Processing Protocol

The simulation output can be processed using a suite of tools in the `postProcessScripts/` directory. These tools generate visualizations and extract quantitative data from the simulation results.

### Quick Start

1. Navigate to your test case directory:
```bash
cd testCases/
```

2. Run the post-processing script:<br>
-> The folderToProcess is the name of the folder containing the simulation results (this folder will contain the intermediate folder with all the simulation snapshots).
```bash
./postProcessData.sh <folderToProcess>
```

### Components and Tools

#### 1. Post-Processing Script (`postProcessData.sh`)
- Main script that orchestrates the post-processing workflow
- Copies required executables and scripts to the working directory
- Executes visualization and data extraction routines
- Cleans up temporary files after processing

#### 2. Visualization (`video.py`)
- Python script for generating visualizations and videos
- Features:
  - Custom color maps for better visualization
  - LaTeX-rendered labels and annotations
  - Automated video generation of wave evolution
- Dependencies:
  - NumPy
  - Matplotlib
  - Subprocess
  - Multiprocessing for parallel processing

#### 3. Data Extraction Tools
- `getData`: Extracts field data (velocity, pressure, etc.)
  - Usage: `./getData filename zmin rmin zmax rmax nr`
  - Outputs: Field values on specified grid points

- `getFacets`: Extracts interface geometry
  - Usage: `./getFacets filename`
  - Outputs: Interface segments as coordinate pairs

### Output Files

The post-processing generates several types of output:
1. Video files showing wave evolution
2. Data files containing:
   - Grid points (X, Y)
   - velocity magnitude field `vel` (saved as vel_rot)
   - vorticity field `omega` (saved as omega_rot)
   - the vof field `f` (saved as f_rot).
   - Field values (`ux`, `uy`) -- saved as (ux_rot, uy_rot)

### Customization

The visualization parameters can be customized in `video.py`:
- Grid resolution
- Color schemes
- Plot dimensions
- Output format and quality

### Troubleshooting

Common issues and solutions:
1. Missing executables: Ensure `getData` and `getFacets` are compiled
2. Python dependencies: Install required packages using pip
3. Permission issues: Make sure all scripts are executable (`chmod +x`)

## Project Team

- Vatsal Sanjay (University of Twente)
  - Email: vatsalsanjay@gmail.com

## Contributing

If you encounter any issues or have suggestions for improvements, please feel free to [open an issue](https://github.com/vatsal-agarwal-20/standing-non-linear-waves/issues) on GitHub. We welcome bug reports, feature requests, and general feedback to help improve the project.

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.

## Citation

If you use this code in your research, please cite:

```bibtex
@software{vatsal_sanjay_2024_14551990,
  author       = {Vatsal Sanjay},
  title        = {comphy-lab/standing-non-linear-waves: Standing Non-Linear Waves – First Official Release -- v1.0},
  year         = {2024},
  publisher    = {Zenodo},
  version      = {v1.0},
  doi          = {10.5281/zenodo.14551990},
  url          = {https://doi.org/10.5281/zenodo.14551990}
}