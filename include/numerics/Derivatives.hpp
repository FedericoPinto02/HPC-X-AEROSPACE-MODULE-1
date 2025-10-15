#include <cmath>
#include <vector>

class Derivatives {
private:
    double dx, dy, dz;
    double nx,ny,nz;
    /**
     * @brief Compute the derivative
     * @param v vector on which derivative is computed
     */
    std::vector<double> derive_x( std::vector<double>& v) const;
public: 

    /**
     * @brief Constructor.
     */
    Derivatives( double dx_, double dy_, double dz_ , double nx_, double ny_, double nz_);

    /**
     * @brief Compute gradient of an input
     * @param p element on which gradient is computed
     */
    std::vector<std::vector<double>> gradient( std::vector<double>& p) const;

    /**
     * @brief Compute divergence of an input
     * @param v vector on which divergence is computed
     */
    std::vector<double> div( std::vector<std::vector<double>>& v) const;

    /**
     * @brief Compute second derivative over a line
     * @param v vector on which divergence is computed
     * @param ds step
     */
    std::vector<double> dd( std::vector<double>& v, double ds) const;
};