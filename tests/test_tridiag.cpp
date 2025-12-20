#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "core/TridiagMat.hpp"

// =============================================================================
// Constructor Tests
// =============================================================================

TEST(TridiagMatTest, Constructor_invalidArgument) {
    // Constructor throws if n < 2 (assuming implementation enforces this)
    EXPECT_THROW(TridiagMat matrix(1), std::invalid_argument);
}

TEST(TridiagMatTest, Constructor_correctCreation) {
    TridiagMat matrix(5);

    // Implementation should resize all vectors to size N
    EXPECT_EQ(matrix.getDiag(0).size(), 5);
    EXPECT_EQ(matrix.getDiag(-1).size(), 5);
    EXPECT_EQ(matrix.getDiag(1).size(), 5);
}

// =============================================================================
// GetSize Tests
// =============================================================================

TEST(TridiagMatTest, GetSize) {
    TridiagMat matrix(6);
    EXPECT_EQ(matrix.getSize(), 6);
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

// =============================================================================
// GetFirstElementFromDiag Tests
// =============================================================================

TEST(TridiagMatTest, GetFirstElementFromDiag_Main) {
    TridiagMat matrix(4);
    matrix.getDiag(0) = {10.0, 20.0, 30.0, 40.0};
    EXPECT_EQ(matrix.getFirstElementFromDiag(0), 10.0);
}

TEST(TridiagMatTest, GetFirstElementFromDiag_Sub) {
    TridiagMat matrix(4);
    matrix.getDiag(-1) = {5.0, 6.0, 7.0, 8.0};
    EXPECT_EQ(matrix.getFirstElementFromDiag(-1), 5.0);
}

TEST(TridiagMatTest, GetFirstElementFromDiag_Sup) {
    TridiagMat matrix(4);
    matrix.getDiag(1) = {9.0, 8.0, 7.0, 6.0};
    EXPECT_EQ(matrix.getFirstElementFromDiag(1), 9.0);
}

TEST(TridiagMatTest, GetFirstElementFromDiag_Invalid) {
    TridiagMat matrix(4);
    EXPECT_THROW(matrix.getFirstElementFromDiag(5), std::invalid_argument);
}

// =============================================================================
// GetLastElementFromDiag Tests
// =============================================================================

TEST(TridiagMatTest, GetLastElementFromDiag_Main) {
    TridiagMat matrix(4);
    // Implementation likely returns diag.at(size-1) for main diagonal
    matrix.getDiag(0) = {1.0, 2.0, 3.0, 4.0};
    EXPECT_EQ(matrix.getLastElementFromDiag(0), 4.0);
}

TEST(TridiagMatTest, GetLastElementFromDiag_Sub) {
    TridiagMat matrix(4);
    // Implementation usually returns subdiag.at(size-2) for off-diagonals
    // Indices: 0->8.0, 1->9.0, 2->10.0, 3->11.0
    // size-2 = index 2 => 10.0
    matrix.getDiag(-1) = {8.0, 9.0, 10.0, 11.0};
    EXPECT_EQ(matrix.getLastElementFromDiag(-1), 10.0);
}

TEST(TridiagMatTest, GetLastElementFromDiag_Sup) {
    TridiagMat matrix(4);
    // Implementation usually returns supdiag.at(size-2) for off-diagonals
    // Indices: 0->5.0, 1->6.0, 2->7.0, 3->0.0
    // size-2 = index 2 => 7.0
    matrix.getDiag(1) = {5.0, 6.0, 7.0, 0.0};
    EXPECT_EQ(matrix.getLastElementFromDiag(1), 7.0);
}

TEST(TridiagMatTest, GetLastElementFromDiag_Invalid) {
    TridiagMat matrix(4);
    EXPECT_THROW(matrix.getLastElementFromDiag(-5), std::invalid_argument);
}