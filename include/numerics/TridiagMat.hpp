#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

/**
 * @class TridiagMat
 * @brief Class representing a tridiagonal matrix and providing access to its elements.
 */
class TridiagMat {
private:
    std::vector<double> diag;     ///< Main diagonal elements
    std::vector<double> subdiag;  ///< Sub-diagonal elements
    std::vector<double> supdiag;  ///< Super-diagonal elements

public:
    /**
     * @brief Constructor that initializes the tridiagonal matrix of size n x n.
     * @param n the size of the matrix
     */
    explicit TridiagMat(size_t n = 2);

    /**
     * @brief Get the matrix size.
     */
    [[nodiscard]] inline size_t getSize() const { return diag.size(); }

    /**
     * @brief Resize the matrix to new size n x n.
     * @param n the new size of the matrix
     */
    void resize(size_t n);

    /**
     * @brief Get the whole diagonal, subdiagonal or supdiagonal
     * @param w number indicating sub (-1), diag (0) or sup (1)
     */
    [[nodiscard]] std::vector<double> &getDiag(int w);

    /// @overload
    [[nodiscard]] const std::vector<double> &getDiag(int w) const;
};
