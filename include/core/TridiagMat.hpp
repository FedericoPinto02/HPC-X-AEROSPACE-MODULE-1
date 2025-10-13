#include <vector>
#include <cmath>


class TridiagMat {
  private:
    std::vector<double> diag;
    std::vector<double> subdiag;
    std::vector<double> supdiag;
    const unsigned int size = 0;

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
    void fillMat(std::vector<double> diag, std::vector<double> subdiag, std::vector<double> supdiag);

    /**
     * @brief Get the matrix size
     */
    const unsigned int getSize() const;

    /**
     * @brief Get the (i,j) element of the matrix
     */
    double get(int i, int j) const;
};