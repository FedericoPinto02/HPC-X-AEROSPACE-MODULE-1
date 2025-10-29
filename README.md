# HPC-X-AEROSPACE-MODULE-1

## 🧩 Workflow Draft

### Main Execution Flow
```text
main
 ├─ input.read()
 ├─ initializer.setup()
 └─ timeIntegrator.run()
```

### Time Integration Loop
```text
timeIntegrator
 └─ timeIntegrator.run()
     ├─ for each timeStep
     │   ├─ viscousStep.run()
     │   ├─ pressureStep.run()
     │   └─ log.write()
```

### Solver Computational Steps

#### Viscous Step
```text
viscousStep
 └─ viscousStep.run()
     ├─ viscousStep.computeG()             # parallel
     ├─ viscousStep.computeXi()            # parallel
     └─ for each velocity linsys           # parallel
         ├─ viscousStep.linsys.fillSys()   # keep parallel
         └─ linsys.solve()                 # keep parallel
```

#### Pressure Step
```text
pressureStep
 └─ pressureStep.run()
     ├─ for each pressure linsys           # parallel
     │   ├─ pressureStep.linsys.fillSys()  # keep parallel
     │   └─ linsys.solve()                 # keep parallel
     └─ pressureStep.updatePressure()      # parallel
```

### Utility and Support Modules
```text
linsys
 ├─ objs: triMatrix, unknown, known
 ├─ linsys.solve()
 └─ linsys.fillSys()

derivatives
 ├─ derivatives.laplacian()
 ├─ derivatives.gradient()
 └─ derivatives.divergence()

initializer
 └─ initializer.setup()

input
 └─ input.read()

log
 └─ log.write()
```

---

## 🧱 Program Structure Draft
'''
.
├── benchmarks
│   └── placeholder.txt
├── CMakeLists.txt
├── data
│   └── config.in
├── .DS_Store
├── examples
│   └── placeholder.txt
├── extern
│   └── googletest
├── .gitignore
├── .gitmodules
├── include
│   ├── core
│   │   ├── Fields.hpp
│   │   ├── Mesh.hpp
│   │   └── TridiagMat.hpp
│   ├── io
│   │   ├── inputReader.hpp
│   │   └── logWriter.hpp
│   ├── numerics
│   │   ├── derivatives.hpp
│   │   └── LinearSys.hpp
│   └── simulation
│       ├── initializer.hpp
│       ├── pressureStep.hpp
│       ├── SimulationContext.hpp
│       └── viscousStep.hpp
├── LICENSE
├── README.md
├── results
│   └── results.out
├── src
│   ├── core
│   │   ├── Fields.cpp
│   │   ├── mesh.cpp
│   │   └── TridiagMat.cpp
│   ├── io
│   │   ├── inputReader.cpp
│   │   └── logWriter.cpp
│   ├── main.cpp
│   ├── numerics
│   │   ├── derivatives.cpp
│   │   └── LinearSys.cpp
│   └── simulation
│       ├── initializer.cpp
│       ├── pressureStep.cpp
│       ├── timeIntegrator.cpp
│       └── viscousStep.cpp
├── tests
│   ├── test_derivatives.cpp
│   ├── test_fields.cpp
│   ├── test_linearSys.cpp
│   └── test_tridiag.cpp
├── test.sh
└── .vscode
    └── settings.json
'''