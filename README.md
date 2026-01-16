# NSBSolver  
### **High-Performance Scientific Computing in Aerospace**  
Politecnico di Milano — Academic Year **2025–2026**  
**Group ID:** 1  

---

## Overview

**NSBSolver** is a high-performance parallel **Computational Fluid Dynamics (CFD)** solver for the simulation of **incompressible Stokes–Brinkman equations**.  
It is designed as a computational kernel for **topology optimization** problems in **thermal management** for **Advanced Air Mobility (AAM)** applications.

The solver is written in modern **C++17** and targets **distributed-memory HPC systems**, leveraging **MPI** parallelism, **direction splitting**, and a **Schur-complement-based domain decomposition** strategy.

---

## Contributors

Developed by **Group 1** under the supervision of **Prof. Franco Auteri**:

- Ettore Cirillo  
- Mattia Gotti  
- Giulio Martella  
- Michele Milani  
- Stefano Pedretti  
- Daniele Piano  
- Federico Pinto  

---

## Features

### Physics & Numerical Methods
- **Equations**: Incompressible Stokes–Brinkman
- **Geometry handling**: Fixed-grid approach inspired by the Immersed Boundary Method (IBM) using porosity penalization
- **Spatial discretization**:  
  - Second-order asymmetric **MAC finite-difference** stencil
- **Temporal discretization**:  
  - **Crank–Nicolson** scheme  
  - **Douglas direction splitting** (unconditionally stable)

### Linear Solvers & Algorithms
- Custom **TDMA (Thomas algorithm)** for tridiagonal systems
- **Schur complement method** for parallel domain decomposition
- Explicit decoupling of interior nodes and inter-subdomain interfaces

---

## High-Performance Computing

- **Parallel model**: MPI with Cartesian topology
- **Domain decomposition**: Overlapping subdomains
- **Memory optimization**:
  - Persistent allocation strategy
  - Structure-of-Arrays (SoA) layout for vector fields
- **Communication**:
  - Non-blocking halo exchanges
- **Performance**:
  - Near-linear strong scaling for low core counts  
  - Measured speed: ≈ **0.48 μs per cell-step**  
  - Speedup at 4 cores: **S₄ ≈ 3.3**

---

## Limitations & Future Work

Current limitations include:

- **Convective term**: Not yet included (architecture already supports extension)
- **Boundary conditions**:
  - Supported: Dirichlet (velocity), Neumann (pressure)
  - Missing: Fully periodic BCs
- **Memory layout**:
  - Fields stored as separate objects  
  - Future improvement: single contiguous memory block to reduce TLB misses
- **I/O**:
  - Output currently serialized on master rank  
  - Parallel I/O planned for large-scale runs

---

## Getting Started

### Prerequisites

- **C++ compiler**: C++17 compliant (GCC / Clang)
- **MPI**: OpenMPI or MPICH
- **CMake**: ≥ 3.10

### Build & Run

```bash
# Build the solver
./build.sh

# Configure input parameters
# Edit data/config.json as needed

# Run the simulation
./run.sh

# Run unit tests
./test.sh
```

---

## Project Structure

```text
.
├── data/                       # Runtime configuration
│   ├── config.json             # Grid, time, physics parameters
│   └── configFunctions.hpp     # Analytical BCs and initial conditions
├── include/
│   ├── core/                   # MPI wrappers, grids, distributed fields
│   ├── io/                     # VTK output, JSON parsing, logging
│   ├── numerics/               # FD stencils, TDMA & Schur solvers
│   └── simulation/             # Time stepping and physics operators
├── src/                        # Source files
├── tests/                      # Unit tests
└── output/                     # VTK results
```

---

## Validation

The solver has been validated through:

1. **Method of Manufactured Solutions (MMS)**  
   - Verified second-order accuracy in space and time

2. **Canonical Laminar Flows**
   - Couette–Poiseuille flow
   - Hagen–Poiseuille flow with Brinkman penalization for circular geometries

3. **Scalability Tests**
   - Near-linear strong scaling up to 4 MPI ranks
   - Robust weak scaling on distributed grids

---

## License & Usage

This code was developed for **academic and research purposes** within the  
*High-Performance Scientific Computing in Aerospace* course at Politecnico di Milano.
