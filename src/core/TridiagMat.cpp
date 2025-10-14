#include <cmath>
#include <stdexcept>
#include <vector>

class TridiagMat {
private:
  std::vector<double> diag;
  std::vector<double> subdiag;
  std::vector<double> supdiag;
  const unsigned int size;

public:
  /**
   * @brief Constructor
   * @param n Size of the matrix (n x n)
   */
  TridiagMat(int n) : size(n) {
    if (n < 2)
      throw std::invalid_argument("Matrix size must be at least 2");
    diag.resize(n, 0.0);
    subdiag.resize(n - 1, 0.0);
    supdiag.resize(n - 1, 0.0);
  }

  /**
   * @brief Fill the matrix
   * @param diag matrix diagonal
   * @param subdiag matrix subdiagonal
   * @param supdiag matrix upper diagonal
   */
  void fillMat(std::vector<double> diag_, std::vector<double> subdiag_,
               std::vector<double> supdiag_) {
    // size checks
    if (diag_.size() != size || subdiag_.size() != (size - 1) ||
        supdiag_.size() != (size - 1))
      throw std::invalid_argument(
          "Invalid vector sizes for tridiagonal matrix");

    // move vectors (no copy is done here :D)
    diag = std::move(diag_);
    subdiag = std::move(subdiag_);
    supdiag = std::move(supdiag_);
  }

  /**
   * @brief Get the matrix size
   */
  const unsigned int getSize() const { return size; }

  /**
   * @brief Get the (i,j) element of the matrix
   */
  double getElement(int i, int j) const {
    if (i < 0 || j < 0 || i >= size || j >= size)
      throw std::out_of_range("Index out of range");

    if (i == j)
      return diag[i];
    if (i == j + 1)
      return subdiag[j];
    if (i + 1 == j)
      return supdiag[i];
    return 0.0; // elements outside the tridiago are zero
  }

  /**
   * @brief Get the whole diagonal, subdiagonal or supdiagonal
   * @param w number indicating sub (-1), diag (0) or sup (1)
   */
  std::vector<double> getDiags(int w) const {
    if (w == -1)
      return subdiag;
    else if (w == 0)
      return diag;
    else if (w == 1)
      return supdiag;
    else {
      throw std::invalid_argument("Parameter needs to be -1, 0 or 1");
    }
  }

  /**
   * @brief Get the first element from the diagonal w
   * @param w number indicating sub (-1), main (0) or sup (1) diagonal
   */
  double getFirstElementFromDiag(int w) const {
    if (w == -1)
      return subdiag.at(0);
    else if (w == 0)
      return diag.at(0);
    else if (w == 1)
      return supdiag.at(0);
    else
      throw std::invalid_argument("Parameter needs to be -1, 0 or 1");
  }

  /**
   * @brief Get the last element from the diagonal w
   * @param w number indicating sub (-1), main (0) or sup (1) diagonal
   */
  double getLastElementFromDiag(int w) const {
    if (w == -1)
      return subdiag.at(size - 2);
    else if (w == 0)
      return diag.at(size - 1);
    else if (w == 1)
      return supdiag.at(size - 2);
    else
      throw std::invalid_argument("Parameter needs to be -1, 0 or 1");
  }
};