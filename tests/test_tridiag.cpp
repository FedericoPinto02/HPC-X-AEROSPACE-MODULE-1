#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "numerics/TridiagMat.hpp"

// =============================================================================
// Constructor Tests
// =============================================================================

TEST(TridiagMatTest, Constructor_invalidArgument) {
    // Constructor throws if n < 2 (assuming implementation enforces this)
    EXPECT_THROW(TridiagMat matrix(1), std::invalid_argument);
}

TEST(TridiagMatTest, Constructor_correctCreation) {
    TridiagMat matrix(5);

    size_t subdiag_size = matrix.getDiag(-1).size();
    size_t diag_size    = matrix.getDiag(0).size();
    size_t supdiag_size = matrix.getDiag(1).size();

    // Implementation should resize all vectors to size N
    EXPECT_EQ(diag_size, 5);
    EXPECT_EQ(subdiag_size, 5);
    EXPECT_EQ(supdiag_size, 5);
}

// =============================================================================
// GetSize Tests
// =============================================================================

TEST(TridiagMatTest, GetSize) {
    TridiagMat matrix(6);
    size_t size = matrix.getSize();
    EXPECT_EQ(size, 6);
}

// =============================================================================
// GetDiag Tests (Read/Write)
// =============================================================================

TEST(TridiagMatTest, GetDiag_WriteAndRead) {
    TridiagMat matrix(4);

    // Create dummy data
    std::vector<double> main = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> sub  = {8.0, 9.0, 10.0, 11.0};
    std::vector<double> sup  = {5.0, 6.0, 7.0, 0.0};

    // Write via non-const reference
    matrix.getDiag(0) = main;
    matrix.getDiag(-1) = sub;
    matrix.getDiag(1) = sup;

    // Verify
    EXPECT_EQ(matrix.getDiag(0), main);
    EXPECT_EQ(matrix.getDiag(-1), sub);
    EXPECT_EQ(matrix.getDiag(1), sup);
}

TEST(TridiagMatTest, GetDiag_InvalidArgument) {
    TridiagMat matrix(4);
    // Should throw if parameter is not -1, 0, or 1
    EXPECT_THROW(matrix.getDiag(-2), std::invalid_argument);
    EXPECT_THROW(matrix.getDiag(2), std::invalid_argument);
}

TEST(TridiagMatTest, GetDiag_ModifyReference) {
    TridiagMat matrix(4);
    std::vector<double> sub = {8.0, 9.0, 10.0, 11.0};
    matrix.getDiag(-1) = sub;

    // Modify via reference
    std::vector<double> &ref = matrix.getDiag(-1);
    ref.front() *= 2; // 16.0
    ref.back() *= 2;  // 22.0

    EXPECT_DOUBLE_EQ(matrix.getDiag(-1)[0], 16.0);
    EXPECT_DOUBLE_EQ(matrix.getDiag(-1)[3], 22.0);
}
