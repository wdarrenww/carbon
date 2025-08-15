# Carbon – Molecular Interaction Simulator in C++

Carbon is a C++ simulation engine for constructing and observing molecular systems with a focus on physical realism.  
It is designed as an educational and research-oriented tool to explore how atoms and molecules interact under classical force fields.

---

## Features

- Manual or file-based molecule construction
- Bonded and nonbonded interactions:
  - Harmonic bonds, angles, and dihedrals
  - Lennard-Jones and Coulomb potentials
- Velocity-Verlet integration
- Neighbor list optimization for performance
- Periodic boundary conditions
- Langevin thermostat (optional)
- XYZ trajectory output for visualization
- Extensible force field and analysis modules

---

## Project Goals

1. Provide a modular, extensible simulation core.
2. Balance computational efficiency with physical realism.
3. Allow interactive construction and modification of molecular systems.
4. Serve both educational and research use cases.

---

## Directory Overview

- **`include/`** – Public header files organized by subsystem.
- **`src/`** – Source implementations matching the `include` structure.
- **`data/`** – Example molecules, parameter files, and sample systems.
- **`tests/`** – Unit tests for core physics and I/O.
- **`docs/`** – Technical documentation and design notes.

---

## Build Instructions

### Prerequisites
- C++18 or later
- CMake 3.16+
- Eigen (for vector/matrix math)
- nlohmann/json (for parameter files)
- Optional: OpenGL + Dear ImGui for GUI visualization

---

## Parameter Files

Parameters are stored in JSON format (data/parameters.json) with per-element values

---

## Roadmap
See docs/design.md for a phased development plan:
1. Core engine and basic force terms.
2. Neighbor lists, thermostats, and stability.
3. Angles, dihedrals, and more realistic force fields.
4. Interactive GUI for building and viewing molecules.
5. Optional advanced features (reactive force fields, QM/MM coupling).

---

## License
Licensed under the MIT License. See LICENSE for details.

---

## Contributing
Contributions are welcome.
Please read docs/design.md for code style and architecture guidelines before submitting pull requests.

---