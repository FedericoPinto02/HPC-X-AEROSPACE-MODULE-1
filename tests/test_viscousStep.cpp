#include <gtest/gtest.h>
#include <memory>
#include <vector>

// --- Include headers ---
#include "simulation/viscousStep.hpp"
#include "simulation/SimulationContext.hpp" // Assuming this exists
#include "core/Fields.hpp"
#include "numerics/derivatives.hpp"

/// In tests/test_viscousStep.cpp
// In tests/test_viscousStep.cpp

class ViscousStepLinkTest : public ::testing::Test {
protected:
    // Test objects
    std::shared_ptr<Grid> grid;
    SimulationContext context;
    std::unique_ptr<ViscousStep> viscousStep;

    // --- CONSTRUCTOR ---
    ViscousStepLinkTest()
        : // 1. First, create the grid
          grid(std::make_shared<Grid>(3, 3, 3, 1.0, 1.0, 1.0)),
          // 2. Now, initialize 'context' in the correct order:
          //    - 1st member: gridPtr (use the 'grid' we just made)
          //    - 2nd member: timeSettings (use {dt, totalTime})
          context{ grid, {0.1, 1.0} }
    {
        // Constructor body can be empty.
        // SetUp() will run immediately after this.
    }

    // SetUp() runs after the constructor
    // In test_viscousStep.cpp
    void SetUp() override {
        // ...
        std::vector<Field::Scalar> zeros(grid->size(), 0.0);
        std::vector<Field::Scalar> ones(grid->size(), 1.0); 

        // ...
        context.constants.nu = 1.0;
        context.constants.k.setup(grid, ones); 
        context.constants.f.setup(grid, zeros, zeros, zeros);
        context.state.p.setup(grid, zeros);

        // Inizializza sia i campi 'attuali' SIA quelli 'vecchi'
        context.state.u.setup(grid, zeros, zeros, zeros);
        context.state.eta.setup(grid, zeros, zeros, zeros);
        context.state.zeta.setup(grid, zeros, zeros, zeros);
        
        context.state.uOld.setup(grid, zeros, zeros, zeros);
        context.state.etaOld.setup(grid, zeros, zeros, zeros);
        context.state.zetaOld.setup(grid, zeros, zeros, zeros);

        // 4. Create the object under test
        viscousStep = std::make_unique<ViscousStep>(context);
    }
};

// ... (The TEST_F functions remain unchanged) ...


// --- The Tests ---

// Test 1: Check if constructor and initialization work
TEST_F(ViscousStepLinkTest, ConstructorAndInit) {
    // If SetUp() succeeded, viscousStep is not null.
    ASSERT_NE(viscousStep, nullptr);
}

// Test 2: Check if computeG() compiles, links, and does not crash
TEST_F(ViscousStepLinkTest, ComputeGLinks) {
    // This test will fail at link-time if 
    // computeG OR computeGradient, Dxx, Dyy, Dzz are not defined.
    ASSERT_NO_THROW(viscousStep->run());
}
