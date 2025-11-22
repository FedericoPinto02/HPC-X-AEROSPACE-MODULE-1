# HPC-X-AEROSPACE-MODULE-1: NSBSolver (Navier–Stokes–Brinkman Solver)

NSBSolver is a high-performance Computational Fluid Dynamics (CFD) module designed to solve the **incompressible Navier–Stokes–Brinkman equations**.

The solver utilizes a **Fractional-Step Projection Method** for time integration and features a custom **Domain Decomposition** linear solver based on the **Schur Complement method** to handle large-scale linear systems sequentially.

---

## 🛠️ Quick Start

### Prerequisites
* **C++ Compiler**: C++17 compliant (GCC, Clang).
* **CMake**: Version 3.10 or higher.

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

## 🚀 Key Features

### Core Physics
* **Governing Equations**: Incompressible Navier–Stokes–Brinkman.
* **Domain Types**: Handles transitions between pure fluid (Navier-Stokes) and porous media (Brinkman/Darcy) via variable permeability fields.

### Numerical Method
* **Spatial Discretization**: Finite Difference Method (FDM) on a Staggered Grid.
* **Time Integration**: Fractional-step Projection Method:
    1.  **Viscous Step**: Solves the momentum equation.
    2.  **Pressure Step**: Solves the Poisson equation.

---


## 📂 Project Structure

```text
.
├── data/                       # Configuration
│   ├── config.json             # Runtime parameters
│   └── configFunctions.hpp     # Boundary Condition definitions
├── include/
│   ├── core/                   # Mesh, Fields, TridiagMat
│   ├── io/                     # VTKWriter, InputReader, LogWriter
│   ├── numerics/               # LinearSys, SchurSequentialSolver, Derivatives
│   └── simulation/             # NSBSolver, ViscousStep, PressureStep
├── src/                        # Source implementation
├── tests/                      # Unit testing suite
└── output/                     # VTK simulation results
```

---

## 🧩 Workflow Architecture

The simulation is orchestrated by the `NSBSolver` class, which manages the lifecycle of the simulation data and the time-stepping loop.

### Main Execution Flow

```text
main()
 └─ NSBSolver solver("config.json")
     ├─ setup()
     │   ├─ InputReader::read()          # Parse JSON
     │   ├─ Initializer::setup()         # Allocate Grids, Fields, and BCs
     │
     └─ solve()                          # Main Time Loop
         ├─ ViscousStep::run()           # Predictor
         ├─ PressureStep::run()          # Corrector
         └─ VTKWriter::write()           # Visualization export
```