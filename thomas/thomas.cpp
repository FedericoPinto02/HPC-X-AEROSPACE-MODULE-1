#include <vector>
#include <stdexcept>
using namespace std;

void thomas( vector<double>& a, vector<double>& b , vector<double>& c, 
                vector<double>& f, vector<double>& x )
{
    /** 
    * @brief Solves a tridiagonal system Ax = f using the Thomas algorithm.
    *
    * @param a Sub-diagonal coefficients (size n-1)
    * @param b Main diagonal coefficients (size n)
    * @param c Super-diagonal coefficients (size n-1)
    * @param f Right-hand side vector (size n)
    * @param x Output vector containing the solution (size n)
    *
    * @throws std::invalid_argument if the vector sizes are not consistent.
    */

    // Compatibility check
    const unsigned int n = b.size();
    if (a.size() != n - 1 || c.size() != n - 1 || f.size() != n)
        throw invalid_argument("Dimension mismatch for a tridiagonal system");

    // Forward Sweep
    c[0] /= b[0];
    f[0] /= b[0];
    for ( unsigned int i = 1; i < n ; i++)
    {
        b[i] = b[i] - a[i-1] * c[i-1];
        f[i] = f[i] - a[i-1] * f[i-1];
        if (i < n - 1)                  // needed because of c dimensions (n-1)
            c[i] /= b[i];
        f[i] /= b[i];
    }

    // Backward Sweep 
    x.resize(n);
    x[n-1] = f[n-1];
    for ( int i = n - 2; i >= 0; i-- )
    {
        x[i] = f[i] - c[i] * x[i+1];
    }
}

