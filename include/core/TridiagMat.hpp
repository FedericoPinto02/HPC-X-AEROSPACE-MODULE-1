#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

///Class representing a tridiagonal matrix and providing access to its elements.
class TridiagMat {
private:
    std::vector<double> diag;
    std::vector<double> subdiag;
    std::vector<double> supdiag;

public:
    /**
     * @brief Constructor.
     * @param n Size of the matrix (n x n)
     */
    explicit TridiagMat(size_t n);

    /**
     * @brief Fill the matrix.
     * @param diag the matrix diagonal
     * @param subdiag the matrix subdiagonal
     * @param supdiag the matrix upper diagonal
     */
    void fillMat(std::vector<double> diag, std::vector<double> subdiag, std::vector<double> supdiag);

    /**
     * @brief Get the matrix size.
     */
    [[nodiscard]] inline size_t getSize() const { return diag.size(); }

    /**
     * @brief Get the (i,j) element of the matrix.
     */
    [[nodiscard]] double getElement(size_t i, size_t j) const;

    /**
     * @brief Get the whole diagonal, subdiagonal or supdiagonal
     * @param w number indicating sub (-1), diag (0) or sup (1)
     */
    [[nodiscard]] std::vector<double> &getDiag(int w);

    /// @overload
    [[nodiscard]] std::vector<double> getDiag(int w) const;

    /**
     * @brief Get the first element from the diagonal w
     * @param w number indicating sub (-1), main (0) or sup (1) diagonal
     */
    [[nodiscard]] double getFirstElementFromDiag(int w) const;

    /**
     * @brief Get the last element from the diagonal w
     * @param w number indicating sub (-1), main (0) or sup (1) diagonal
     */
    [[nodiscard]] double getLastElementFromDiag(int w) const;
};
