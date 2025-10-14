#include <cmath>
#include <vector>

class TridiagMat {
private:
  std::vector<double> diag;
  std::vector<double> subdiag;
  std::vector<double> supdiag;
  const unsigned int size;

public:
  /**
   * @brief Constructor.
   * @param n Size of the matrix (n x n)
   */
  TridiagMat(int n);

  /**
   * @brief Fill the matrix.
   * @param diag matrix diagonal.
   * @param subdiag matrix subdiagonal.
   * @param supdiag matrix upper diagonal.
   */
  void fillMat(std::vector<double> diag, std::vector<double> subdiag,
               std::vector<double> supdiag);

  /**
   * @brief Get the matrix size
   */
  const unsigned int getSize() const;

  /**
   * @brief Get the (i,j) element of the matrix
   */
  double getElement(int i, int j) const;

  /**
   * @brief Get the whole diagonal, subdiagonal or supdiagonal
   * @param w number indicating sub (-1), diag (0) or sup (1)
   */
  std::vector<double> getDiags(int w) const;

  /**
   * @brief Get the first element from the diagonal w
   * @param w number indicating sub (-1), main (0) or sup (1) diagonal
   */
  double getFirstElementFromDiag(int w) const;

  /**
   * @brief Get the last element from the diagonal w
   * @param w number indicating sub (-1), main (0) or sup (1) diagonal
   */
  double getLastElementFromDiag(int w) const;
};