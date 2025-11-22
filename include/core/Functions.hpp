#ifndef NSBSOLVER_FUNCTIONS_HPP
#define NSBSOLVER_FUNCTIONS_HPP

#include <functional>

namespace Functions {
    using Func = std::function<double(double x, double y, double z, double t)>;

    const Func ZERO = [](double /*x*/, double /*y*/, double /*z*/, double /*t*/ = 0) {
        return 0.0;
    };
} // namespace Functions

#endif //NSBSOLVER_FUNCTIONS_HPP