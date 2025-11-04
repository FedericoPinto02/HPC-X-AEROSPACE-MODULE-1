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
> Files marked with `$$` are currently empty placeholders.

```text
.
├── build
│   └── cmake_file $$
├── CMakeLists.txt $$
├── data
│   └── config.in $$
├── include
│   ├── core
│   │   ├── mesh.hpp $$
│   │   ├── physicalField.hpp
│   │   └── TridiagMat.hpp
│   ├── io
│   │   ├── inputReader.hpp $$
│   │   └── logWriter.hpp $$
│   ├── numerics
│   │   ├── derivatives.hpp $$
│   │   └── linearsys.hpp $$
│   └── simulation
│       ├── initializer.hpp $$
│       ├── pressureStep.hpp $$
│       ├── timeIntegrator.hpp $$
│       └── viscousStep.hpp $$
├── LICENSE
├── README.md
├── results
│   └── results.out $$
└── src
    ├── core
    │   ├── mesh.cpp $$
    │   ├── physicalField.cpp $$ → write .cpp
    │   └── TridiagMat.cpp $$ → write .cpp
    ├── io
    │   ├── inputReader.cpp $$
    │   └── logWriter.cpp $$
    ├── main.cpp $$
    ├── numerics
    │   ├── derivatives.cpp $$
    │   ├── linearsys.cpp $$
    │   ├── thomas.cpp → to merge into linearsys.cpp
    │   └── thomas_test.cpp → to merge into linearsys.cpp
    └── simulation
        ├── initializer.cpp $$
        ├── pressureStep.cpp $$
        ├── timeIntegrator.cpp $$
        └── viscousStep.cpp $$
└── tests
```
