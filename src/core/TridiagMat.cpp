#include <vector>
#include <cmath>

class TridiagMat{

private:

std::vector<double> diag;
std::vector<double> subdiag;
std::vector<double> supdiag;
int size = 0;

public:

/**
  * @brief Constructor.
 */
TridiagMat(int n);



/**
 * @brief Fill the matrix.
 * @param diag matrix diagonal.
 * @param subdiag matrix subdiagonal.
 * @param supdiag matrix upper diagonal.
 */
void fillMat(std::vector<double> diag, std::vector<double> subdiag, std::vector<double> supdiag);


};