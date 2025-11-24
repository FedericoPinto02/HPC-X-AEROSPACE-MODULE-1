#include <cmath>
#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <stdexcept> // Required for std::runtime_error

// Include necessary class definitions
#include "core/Grid.hpp"     // For Grid
#include "core/Fields.hpp"   // For Field and VectorField
#include "numerics/LinearSys.hpp" // The class we are testing

/**
 * @brief Standalone test fixture for the LinearSys class.
 *
 * This fixture does not inherit from any other test fixture.
 * It creates all necessary objects (Grid, Field, VectorField)
 * from scratch required to test the LinearSys interface.
 *
 * It is declared as a 'friend' of LinearSys to allow access
 * to private members for testing.
 */
class LinearSysTestFixture : public ::testing::Test {
protected:
    // --- Grid Configuration ---
    // Using sizes greater than 2 for robust testing
    const size_t Nx = 3, Ny = 4, Nz = 5;
    const double dx = 0.1, dy = 0.1, dz = 0.1;
    const size_t vectorSize = Nx * Ny * Nz; // 60
    std::shared_ptr<const Grid> gridPtr;

    // --- Fields required for testing ---
    Field testField;    // Used for fillSystemPressure
    Field porosity;     // Used for fillSystemVelocity
    VectorField eta, xi; // Used for fillSystemVelocity

    /**
     * @brief Sets up all required objects before each test.
     */
    void SetUp() override {
        // 1. Initialize the Grid
        gridPtr = std::make_shared<const Grid>(Nx, Ny, Nz, dx, dy, dz);

        // 2. Initialize scalar Fields (testField and porosity)
        std::vector<Field::Scalar> scalarData(vectorSize, 1.0);
        
        std::vector<Field::Scalar> testFieldData = scalarData;
        testField.setup(gridPtr, testFieldData);

        std::vector<Field::Scalar> porosityData = scalarData;
        porosity.setup(gridPtr, porosityData);

        // 3. Initialize VectorFields (eta, xi, uBoundNew, uBoundOld)
        std::vector<Field::Scalar> vecDataX(vectorSize, 1.1);
        std::vector<Field::Scalar> vecDataY(vectorSize, 2.2);
        std::vector<Field::Scalar> vecDataZ(vectorSize, 3.3);

        std::vector<Field::Scalar> eta_x = vecDataX, eta_y = vecDataY, eta_z = vecDataZ;
        eta.setup(gridPtr, eta_x, eta_y, eta_z);

        std::vector<Field::Scalar> xi_x = vecDataX, xi_y = vecDataY, xi_z = vecDataZ;
        xi.setup(gridPtr, xi_x, xi_y, xi_z);
    }

    /**
     * @brief Helper function to set up the 3x3 test matrix.
     * This function is a member of LinearSysTestFixture, which IS a friend
     * of LinearSys and can access its private members.
     */
    void setupTestMatrix(LinearSys& sys) {
        // Access private matA
        sys.matA.getDiag( 0).at(0) = 2.0;
        sys.matA.getDiag( 0).at(1) = 2.0;
        sys.matA.getDiag( 0).at(2) = 2.0;
        
        sys.matA.getDiag(-1).at(0) = -1.0;
        sys.matA.getDiag(-1).at(1) = -1.0;
        
        sys.matA.getDiag( 1).at(0) = -1.0;
        sys.matA.getDiag( 1).at(1) = -1.0;
    }
    
    // No TearDown() is needed
};


// === CONSTRUCTOR TESTS ===

TEST_F(LinearSysTestFixture, Constructor_CanBeCreated) {
    // Test that the constructor doesn't throw exceptions with valid sizes
    EXPECT_NO_THROW(LinearSys sys(Nx, BoundaryType::Normal));
    EXPECT_NO_THROW(LinearSys sys(Ny, BoundaryType::Tangent));
    EXPECT_NO_THROW(LinearSys sys(Nz, BoundaryType::Normal));
}

// === SETTER/GETTER TESTS ===

TEST_F(LinearSysTestFixture, SetRhs_MismatchedRhsSize) {
    // System is declared with size Nx
    LinearSys linearSys(Nx, BoundaryType::Normal);
    
    // But the rhs vector has the WRONG size (Nx + 1)
    std::vector<double> rhs(Nx + 1); 

    // Expect a runtime_error because rhs size (Nx+1) doesn't match
    // the system size (Nx)
    EXPECT_THROW(
        linearSys.setRhs(rhs),
        std::runtime_error);
}

// === fillSystemPressure TESTS ===

TEST_F(LinearSysTestFixture, FillSystemPressure_CorrectExecution) {


    // Test for X-direction
    LinearSys linearSysX(Nx, BoundaryType::Normal);
    // Set the known RHS vector using the public setter
    std::vector<double> rhsx = {1.0, 0.0, 1.0};
    linearSysX.setRhs(rhsx);
    // The function now fills the internal rhsC, no parameter needed
    EXPECT_NO_THROW(
        linearSysX.fillSystemPressure(testField, Axis::X));

    // Test for Y-direction
    LinearSys linearSysY(Ny, BoundaryType::Normal);
    std::vector<double> rhsy = {1.0, 0.0, 1.0, 2.0};
    linearSysY.setRhs(rhsy);
    EXPECT_NO_THROW(
        linearSysY.fillSystemPressure(testField, Axis::Y));
    
    // Test for Z-direction
    LinearSys linearSysZ(Nz, BoundaryType::Normal);
    std::vector<double> rhsz = {1.0, 0.0, 1.0, 2.0, 0.5};
    linearSysZ.setRhs(rhsz);
    EXPECT_NO_THROW(
        linearSysZ.fillSystemPressure(testField, Axis::Z));
}


// === fillSystemVelocity TESTS ===

TEST_F(LinearSysTestFixture, FillSystemVelocity_CorrectExecution) {


    // System declared with correct size (Nz) for Axis::Z derivative
    LinearSys linearSys(Nz, BoundaryType::Tangent);
     std::vector<double> rhsz = {1.0, 0.0, 1.0, 2.0, 0.5};
    linearSys.setRhs(rhsz);
    double nu = 1.0e-5;
    double dt = 1.0e-3;


    // The function now fills the internal rhsC, no parameter needed
    EXPECT_NO_THROW(
        linearSys.fillSystemVelocity(porosity, eta, xi, uBoundNew, uBoundOld,
                                     Axis::X, Axis::Z, 0, 0, 0, nu, dt));
    
    // System declared with correct size (Ny) for Axis::Y derivative
    LinearSys linearSysY(Ny, BoundaryType::Tangent);
     std::vector<double> rhsy = {1.0, 0.0, 1.0, 2.0};
    linearSysY.setRhs(rhsy);
    
    EXPECT_NO_THROW(
        linearSysY.fillSystemVelocity(porosity, eta, xi, uBoundNew, uBoundOld,
                                     Axis::Y, Axis::Y, 0, 0, 0, nu, dt));
}

// === SOLVER TESTS ===

/**
 * @brief Tests the ThomaSolver with a manually constructed system.
 */
TEST_F(LinearSysTestFixture, ThomaSolver_CalculatesKnownSolution) {
    const int n = 3;
    LinearSys sys(n, BoundaryType::Normal);

    // 1. Manually create a simple 3x3 system
    //    [ 2 -1  0 ] [x0]   [ 1]
    //    [-1  2 -1 ] [x1] = [ 0]
    //    [ 0 -1  2 ] [x2]   [ 1]
    //    Solution is: x = [1, 1, 1]
    
    // Call the helper function (which has friend access)
    setupTestMatrix(sys); 

    // 2. Set the known RHS vector using the public setter
    std::vector<double> rhs = {1.0, 0.0, 1.0};
    sys.setRhs(rhs);

    // 3. Solve the system
    sys.ThomaSolver();

    // 4. Get the solution vector using the public getter
    const std::vector<double>& solution = sys.getSolution();

    // 5. Check the results
    ASSERT_EQ(solution.size(), n);
    EXPECT_NEAR(solution.at(0), 1.0, 1e-9);
    EXPECT_NEAR(solution.at(1), 1.0, 1e-9);
    EXPECT_NEAR(solution.at(2), 1.0, 1e-9);
}


// === ADDITIONAL TESTS ===

/**
 * @brief Tests the constructor's handling of invalid (zero or negative) sizes.
 *
 * The LinearSys class should not be constructible with a size
 * less than or equal to zero. It should throw an exception.
 */
TEST_F(LinearSysTestFixture, Constructor_ThrowsOnInvalidSize) {
    // A system size of 0 is invalid
    EXPECT_THROW(
        LinearSys sys(0, BoundaryType::Normal),
        std::exception); // Expect any std::exception (or be more specific if known)

    // A negative system size is invalid
    EXPECT_THROW(
        LinearSys sys(-5, BoundaryType::Normal),
        std::exception);
}

/**
 * @brief Tests the getSolution() method on a newly constructed system.
 *
 * When a LinearSys object is first created, the internal solution
 * vector (unknownX) should exist, have the correct size 'n',
 * and be initialized (e.g., to all zeros).
 */
TEST_F(LinearSysTestFixture, GetSolution_ReturnsCorrectSizedVector) {
    const int n = 5;
    LinearSys sys(n, BoundaryType::Normal);
    
    // Get the solution vector before solving
    const std::vector<double>& solution = sys.getSolution();

    // Check that the vector has the correct size
    ASSERT_EQ(solution.size(), n);

    // Check that all initial values are 0.0
    // (This assumes default std::vector<double>(n) initialization)
    for (double val : solution) {
        EXPECT_DOUBLE_EQ(val, 0.0);
    }
}


/**
 * @brief Tests the ThomaSolver with the same matrix but a different RHS.
 *
 * This test uses the same 3x3 matrix from 'setupTestMatrix' but
 * provides a different RHS vector to produce a different
 * known solution.
 */
TEST_F(LinearSysTestFixture, ThomaSolver_CalculatesDifferentSolution) {
    const int n = 3;
    LinearSys sys(n, BoundaryType::Normal);

    // 1. Use the same 3x3 matrix:
    //    [ 2 -1  0 ]
    //    [-1  2 -1 ]
    //    [ 0 -1  2 ]
    setupTestMatrix(sys);

    // 2. Set a new RHS vector: b = {0, 0, 4}
    //    For the matrix A, A*x = b with x = {1, 2, 3} gives:
    //    b_0 = 2(1) - 1(2) + 0(3) = 0
    //    b_1 = -1(1) + 2(2) - 1(3) = -1 + 4 - 3 = 0
    //    b_2 = 0(1) - 1(2) + 2(3) = -2 + 6 = 4
    std::vector<double> rhs = {0.0, 0.0, 4.0};
    sys.setRhs(rhs);

    // 3. Solve the system
    sys.ThomaSolver();

    // 4. Get the solution
    const std::vector<double>& solution = sys.getSolution();

    // 5. Check the new solution: x = {1, 2, 3}
    ASSERT_EQ(solution.size(), n);
    EXPECT_NEAR(solution.at(0), 1.0, 1e-9);
    EXPECT_NEAR(solution.at(1), 2.0, 1e-9);
    EXPECT_NEAR(solution.at(2), 3.0, 1e-9);
}

/**
 * @brief Tests setRhs with another mismatched size (empty vector).
 *
 * This complements the existing 'SetRhs_MismatchedRhsSize' test
 * by checking the edge case of an empty vector.
 */
TEST_F(LinearSysTestFixture, SetRhs_MismatchedRhsSize_EmptyVector) {
    // System is declared with size Nx (3)
    LinearSys linearSys(Nx, BoundaryType::Normal);
    
    // Pass an empty vector
    std::vector<double> rhs_empty; 

    // Expect a runtime_error because rhs size (0) doesn't match
    // the system size (Nx)
    EXPECT_THROW(
        linearSys.setRhs(rhs_empty),
        std::runtime_error);
}