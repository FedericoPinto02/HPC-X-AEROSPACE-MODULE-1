#include <gtest/gtest.h>

#include <cmath>

#include "core/TridiagMat.hpp"

// === Constructor tests ===

TEST(TridiagMatTest, Constructor_invalidArgument) {
    EXPECT_THROW(TridiagMat matrix{1}, std::invalid_argument);
}

TEST(TridiagMatTest, Constructor_correctCreation) {
    TridiagMat matrix(5);

    EXPECT_EQ(matrix.getDiag(0).size(), 5);
    EXPECT_EQ(matrix.getDiag(-1).size(), 4);
    EXPECT_EQ(matrix.getDiag(1).size(), 4);
}

// === fillMat tests ===

TEST(TridiagMatTest, FillMat_incorrectSizeDiag) {
    TridiagMat matrix(3);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0};
    std::vector<double> subdiag = {7.0, 8.0};

    EXPECT_THROW(matrix.fillMat(diag, subdiag, supdiag), std::invalid_argument);
}

TEST(TridiagMatTest, FillMat_incorrectSizeSubDiag) {
    TridiagMat matrix(3);
    std::vector<double> diag = {1.0, 2.0, 3.0};
    std::vector<double> supdiag = {5.0, 6.0};
    std::vector<double> subdiag = {7.0};

    EXPECT_THROW(matrix.fillMat(diag, subdiag, supdiag), std::invalid_argument);
}

TEST(TridiagMatTest, FillMat_incorrectSizeSupDiag) {
    TridiagMat matrix(3);
    std::vector<double> diag = {1.0, 2.0, 3.0};
    std::vector<double> supdiag = {5.0, 6.0, 6.5};
    std::vector<double> subdiag = {7.0, 8.0};

    EXPECT_THROW(matrix.fillMat(diag, subdiag, supdiag), std::invalid_argument);
}

TEST(TridiagMatTest, FillMat_correctPopulation) {
    TridiagMat matrix(3);
    std::vector<double> diag = {1.0, 2.0, 3.0};
    std::vector<double> supdiag = {4.0, 5.0};
    std::vector<double> subdiag = {6.0, 7.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_EQ(matrix.getDiag(0), std::vector<double>({1.0, 2.0, 3.0}));
    EXPECT_EQ(matrix.getDiag(1), std::vector<double>({4.0, 5.0}));
    EXPECT_EQ(matrix.getDiag(-1), std::vector<double>({6.0, 7.0}));
}

// === getSize test ===

TEST(TridiagMatTest, GetSize) {
    TridiagMat matrix(6);

    EXPECT_EQ(matrix.getSize(), 6);
}

// === getElement tests ===

TEST(TridiagMatTest, GetElement_rowIndexTooSmall) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_THROW(matrix.getElement(-2, 1), std::out_of_range);
}

TEST(TridiagMatTest, GetElement_columnIndexTooSmall) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_THROW(matrix.getElement(2, -1), std::out_of_range);
}

TEST(TridiagMatTest, GetElement_rowIndexTooLarge) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_THROW(matrix.getElement(4, 2), std::out_of_range);
}

TEST(TridiagMatTest, GetElement_columnIndexTooLarge) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_THROW(matrix.getElement(3, 5), std::out_of_range);
}

TEST(TridiagMatTest, GetElement_subDiagElement) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_DOUBLE_EQ(matrix.getElement(1, 0), 8.0);
    EXPECT_DOUBLE_EQ(matrix.getElement(2, 1), 9.0);
    EXPECT_DOUBLE_EQ(matrix.getElement(3, 2), 10.0);
}

TEST(TridiagMatTest, GetElement_mainDiagElement) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_DOUBLE_EQ(matrix.getElement(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(matrix.getElement(1, 1), 2.0);
    EXPECT_DOUBLE_EQ(matrix.getElement(2, 2), 3.0);
    EXPECT_DOUBLE_EQ(matrix.getElement(3, 3), 4.0);
}

TEST(TridiagMatTest, GetElement_supDiagElement) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_DOUBLE_EQ(matrix.getElement(0, 1), 5.0);
    EXPECT_DOUBLE_EQ(matrix.getElement(1, 2), 6.0);
    EXPECT_DOUBLE_EQ(matrix.getElement(2, 3), 7.0);
}

// === getDiag tests ===

TEST(TridiagMatTest, GetDiag_subDiag) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_EQ(matrix.getDiag(-1), std::vector({8.0, 9.0, 10.0}));
}

TEST(TridiagMatTest, GetDiag_mainDiag) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_EQ(matrix.getDiag(0), std::vector({1.0, 2.0, 3.0, 4.0}));
}

TEST(TridiagMatTest, GetDiag_supDiag) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_EQ(matrix.getDiag(1), std::vector({5.0, 6.0, 7.0}));
}

TEST(TridiagMatTest, GetDiag_invalidArgument) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_THROW(matrix.getDiag(-2), std::invalid_argument);
}

TEST(TridiagMatTest, GetDiag_changingValuesSubDiag) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    std::vector<double> &matSubDiag = matrix.getDiag(-1);

    matrix.fillMat(diag, subdiag, supdiag);
    matSubDiag.front() *= 2;
    matSubDiag.back() *= 2;

    EXPECT_DOUBLE_EQ(matrix.getFirstElementFromDiag(-1), 16.0);
    EXPECT_DOUBLE_EQ(matrix.getDiag(-1).at(1), 9.0);
    EXPECT_DOUBLE_EQ(matrix.getLastElementFromDiag(-1), 20.0);
}

TEST(TridiagMatTest, GetDiag_changingValuesMainDiag) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    std::vector<double> &matMainDiag = matrix.getDiag(0);

    matrix.fillMat(diag, subdiag, supdiag);
    matMainDiag.front() *= 2;
    matMainDiag.back() *= 2;

    EXPECT_DOUBLE_EQ(matrix.getFirstElementFromDiag(0), 2.0);
    EXPECT_DOUBLE_EQ(matrix.getDiag(0).at(1), 2.0);
    EXPECT_DOUBLE_EQ(matrix.getDiag(0).at(2), 3.0);
    EXPECT_DOUBLE_EQ(matrix.getLastElementFromDiag(0), 8.0);
}

TEST(TridiagMatTest, GetDiag_changingValuesSupDiag) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    std::vector<double> &matSupDiag = matrix.getDiag(1);

    matrix.fillMat(diag, subdiag, supdiag);
    matSupDiag.front() *= 2;
    matSupDiag.back() *= 2;

    EXPECT_DOUBLE_EQ(matrix.getFirstElementFromDiag(1), 10.0);
    EXPECT_DOUBLE_EQ(matrix.getDiag(1).at(1), 6.0);
    EXPECT_DOUBLE_EQ(matrix.getLastElementFromDiag(1), 14.0);
}

// === getFirstElementFromDiag tests ===

TEST(TridiagMatTest, GetFirstElementFromDiag_subDiag) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_EQ(matrix.getFirstElementFromDiag(-1), 8.0);
}

TEST(TridiagMatTest, GetFirstElementFromDiag_mainDiag) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_EQ(matrix.getFirstElementFromDiag(0), 1.0);
}

TEST(TridiagMatTest, GetFirstElementFromDiag_supDiag) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_EQ(matrix.getFirstElementFromDiag(1), 5.0);
}

TEST(TridiagMatTest, GetFirstElementFromDiag_invalidArgument) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_THROW(matrix.getFirstElementFromDiag(3), std::invalid_argument);
}

// === GetLastElementFromDiag tests ===

TEST(TridiagMatTest, GetLastElementFromDiag_subDiag) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_EQ(matrix.getLastElementFromDiag(-1), 10.0);
}

TEST(TridiagMatTest, GetLastElementFromDiag_mainDiag) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_EQ(matrix.getLastElementFromDiag(0), 4.0);
}

TEST(TridiagMatTest, GetLastElementFromDiag_supDiag) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_EQ(matrix.getLastElementFromDiag(1), 7.0);
}

TEST(TridiagMatTest, GetLastElementFromDiag_invalidArgument) {
    TridiagMat matrix(4);
    std::vector<double> diag = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> supdiag = {5.0, 6.0, 7.0};
    std::vector<double> subdiag = {8.0, 9.0, 10.0};

    matrix.fillMat(diag, subdiag, supdiag);
    EXPECT_THROW(matrix.getLastElementFromDiag(-5), std::invalid_argument);
}