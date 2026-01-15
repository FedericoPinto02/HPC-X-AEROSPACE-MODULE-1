
```markdown
# HPC-X-AEROSPACE-MODULE-1: NSBSolver (Navier–Stokes–Brinkman Solver)

**High-Performance Scientific Computing in Aerospace – Politecnico di Milano** **Group ID:** 1  
**Academic Year:** 2025-2026

NSBSolver is a high-performance parallel Computational Fluid Dynamics (CFD) solver designed to simulate **incompressible Stokes–Brinkman equations**. It is specifically engineered to serve as a computational kernel for topology optimization in thermal management applications for Advanced Air Mobility.

The solver features a highly optimized C++ architecture employing a **Direction Splitting** algorithm and a custom **Domain Decomposition** strategy based on the **Schur Complement method**, enabling efficient execution on distributed memory clusters via MPI.

---

## 👥 Contributors

This project was developed by Group 1 under the supervision of **Prof. Franco Auteri**:

* **Ettore Cirillo**
* **Mattia Gotti**
* **Giulio Martella**
* **Michele Milani**
* **Stefano Pedretti**
* **Daniele Piano**
* **Federico Pinto**

---

## 🚀 Key Features & Capabilities

### Core Physics & Numerics
* **Governing Equations**: Incompressible Stokes–Brinkman equations.
* **Methodology**: Fixed-grid strategy inspired by the Immersed Boundary Method (IBM), utilizing porosity penalization to define solid boundaries.
* **Discretization**:
    * **Space**: Second-order asymmetric Marker-And-Cell (MAC) finite difference stencil.
    * **Time**: Unconditionally stable Crank-Nicolson scheme with Direction Splitting (Douglas scheme).
* **Linear Algebra**:
    * Custom **Thomas Algorithm (TDMA)** solver for tridiagonal systems.
    * **Schur Complement Method** for parallel domain decomposition, decoupling internal nodes from shared interfaces.

### High-Performance Computing (HPC)
* **Parallelism**: Distributed memory parallelization using **MPI** with Cartesian topology and overlapping domain decomposition.
* **Optimization**:
    * Persistent memory management to minimize allocation overhead.
    * Structure of Arrays (SoA) layout for vector fields to maximize cache locality.
    * Non-blocking communication for halo exchanges.
* **Performance**: Validated strictly linear scaling for low core counts ($S_4 \approx 3.3$) and computational efficiency of $\approx 0.48 \, \mu s$ per cell-step.

---

## ⚠️ Limitations & Future Work

While the solver meets strict performance targets, current limitations include:
* **Convective Term**: The non-linear convective term is currently neglected to simplify the implementation, though the architecture supports its inclusion.
* **Boundary Conditions**: Currently supports Dirichlet (velocity) and Neumann (pressure). Periodic boundary conditions are not yet fully implemented.
* **Memory Layout**: Fields are currently stored as separate objects. Future optimization could involve aggregating all fields into a single contiguous memory block to reduce TLB misses.
* **IO Overhead**: Parallel I/O is currently serialized (master rank writes), which may become a bottleneck at very large scales.

---

## 🛠️ Quick Start

### Prerequisites
* **C++ Compiler**: C++17 compliant (GCC, Clang) with OpenMPI support.
* **CMake**: Version 3.10 or higher.
* **MPI**: A working MPI implementation (e.g., OpenMPI, MPICH).

### Building and Running
Use the provided helper scripts to compile and run the simulation:

```bash
# 1. Build the project
./build.sh

# 2. Run the simulation
# Ensure data/config.json is configured properly before running
./run.sh

# 3. Run Tests
# Executes unit tests for derivatives, linear solvers, and physics steps
./test.sh

```

---

## 📂 Project Structure

The codebase follows a modular architecture separating physics, numerics, and system management:

```text
.
├── data/                       # Configuration
│   ├── config.json             # Runtime parameters (grid, time, physics)
│   └── configFunctions.hpp     # Compile-time analytical functions (BCs, Initial Conditions)
├── include/
│   ├── core/                   # MPI Wrappers, Grid topology, Distributed Fields
│   ├── io/                     # VTK Visualization, JSON Parsing, Logging
│   ├── numerics/               # Finite Differences, Tridiagonal & Schur Solvers
│   └── simulation/             # Simulation Data, Time Stepper, Physics Steps
├── src/                        # Source implementation
├── tests/                      # Unit testing suite
└── output/                     # VTK simulation results

```

---

## 📊 Validation & Results

The solver has been rigorously validated:

1. **MMS (Method of Manufactured Solutions)**: Confirmed second-order accuracy in space () and time ().
2. **Laminar Flow Benchmarks**:
* **Couette-Poiseuille Flow**: Validated against analytical profiles.
* **Hagen-Poiseuille (Pipe) Flow**: Tested Brinkman penalization for circular geometries.


3. **Scalability**: Demonstrated near-linear strong scalability up to 4 cores and robust weak scaling on distributed grids.

```

```