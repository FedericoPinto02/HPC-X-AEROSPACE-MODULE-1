#include "core/TridiagMat.hpp"
#include <cmath>
#include <stdexcept>
#include <vector>

TridiagMat::TridiagMat(int n) : size(n) {
  if (n < 2)
    throw std::invalid_argument("Matrix size must be at least 2");
  diag.resize(n, 0.0);
  subdiag.resize(n - 1, 0.0);
  supdiag.resize(n - 1, 0.0);
}

void TridiagMat::fillMat(std::vector<double> diag_,
                         std::vector<double> subdiag_,
                         std::vector<double> supdiag_) {
  // size checks
  if (diag_.size() != size || subdiag_.size() != (size - 1) ||
      supdiag_.size() != (size - 1))
    throw std::invalid_argument("Invalid vector sizes for tridiagonal matrix");

  // move vectors (no copy is done here :D)
  diag = std::move(diag_);
  subdiag = std::move(subdiag_);
  supdiag = std::move(supdiag_);
}


double TridiagMat::getElement(int i, int j) const {
  if (i < 0 || j < 0 || i >= size || j >= size)
    throw std::out_of_range("Index out of range");

  if (i == j)
    return diag[i];
  if (i == j + 1)
    return subdiag[j];
  if (i + 1 == j)
    return supdiag[i];
  return 0.0; // elements outside the tridiagonals are zero
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
    return subdiag.at(size - 2);
  else if (w == 0)
    return diag.at(size - 1);
  else if (w == 1)
    return supdiag.at(size - 2);
  else
    throw std::invalid_argument("Parameter needs to be -1, 0 or 1");
}
