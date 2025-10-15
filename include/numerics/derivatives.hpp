#include <cmath>
#include <vector>
#include <string>

class Derivatives {
    private:
        /**
         * @brief Compute the derivative
         * @param v vector on which derivative is computed
         */
        std::vector<double> derive( std::vector<double> v ) const;
    public: 
        /**
         * @brief Constructor.
         */
        Derivatives();

        /**
         * @brief Compute gradient of an input
         * @param p element on which gradient is computed
         */
        std::vector<std::vector<double>> gradient( std::vector<double> p) const;

        /**
         * @brief Compute divergence of an input
         * @param v vector on which divergence is computed
         */
        std::vector<double> div( std::vector<std::vector<double>> v) const;

        /**
         * @brief Compute second derivative over a line
         * @param v vector on which divergence is computed
         */
        std::vector<double> dd( std::vector<double> v) const;
};