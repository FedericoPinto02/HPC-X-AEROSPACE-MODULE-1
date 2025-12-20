#include "core/TridiagMat.hpp"

TridiagMat::TridiagMat(size_t n) {
    if (n < 2)
        throw std::invalid_argument("Matrix size must be at least 2");
    diag.resize(n, 0.0);
    subdiag.resize(n, 0.0);
    supdiag.resize(n, 0.0);
}

std::vector<double> TridiagMat::getDiag(int w) const {
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

std::vector<double> &TridiagMat::getDiag(int w) {
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


double TridiagMat::getFirstElementFromDiag(int w) const {
    if (w == -1)
        return subdiag.at(0);
    else if (w == 0)
        return diag.at(0);
    else if (w == 1)
        return supdiag.at(0);
    else
        throw std::invalid_argument("Parameter needs to be -1, 0 or 1");
}

double TridiagMat::getLastElementFromDiag(int w) const {
    if (w == -1)
        return subdiag.at(diag.size() - 2);
    else if (w == 0)
        return diag.at(diag.size() - 1);
    else if (w == 1)
        return supdiag.at(diag.size() - 2);
    else
        throw std::invalid_argument("Parameter needs to be -1, 0 or 1");
}
