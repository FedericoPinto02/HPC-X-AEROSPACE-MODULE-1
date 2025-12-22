#include "core/Fields.hpp"

#include <stdexcept>

// ---------------------------------------------------------------------------------------------------------------------
// Field class methods
// ---------------------------------------------------------------------------------------------------------------------
void Field::populate(double time) {
    const long Nx = (long) gridPtr_->Nx;
    const long Ny = (long) gridPtr_->Ny;
    const long Nz = (long) gridPtr_->Nz;

    for (long k = 0; k < Nz; ++k) {
        const long kOffset = (long) (k + cached_nHalo_) * (long) cached_axisStrides_[2];
        for (long j = 0; j < Ny; ++j) {
            const long jOffset = (long) (j + cached_nHalo_) * (long) cached_axisStrides_[1];
            for (long i = 0; i < Nx; ++i) {
                auto x = gridPtr_->to_x(i, offset_, offsetAxis_);
                auto y = gridPtr_->to_y(j, offset_, offsetAxis_);
                auto z = gridPtr_->to_z(k, offset_, offsetAxis_);
                data_[kOffset + jOffset + (cached_nHalo_ + i)] = populateFunction_(x, y, z, time);
            }
        }
    }
}

void
Field::reset(const Field::Scalar value) {
    std::fill(data_.begin(), data_.end(), value);
}

void Field::add(const Field::Scalar value) {
    for (auto &elem: data_) {
        elem += value;
    }
    // todo - workaround ??: std::execution::par_unseq not available for clang 15 (libc++ 15)
    /*std::for_each(std::execution::par_unseq, data_.begin(), data_.end(),
                  [value](Scalar &elem) { elem += value; });*/
}

void Field::add(const Field &other) {
    if (other.getGrid().Nx != gridPtr_->Nx
        || other.getGrid().Ny != gridPtr_->Ny
        || other.getGrid().Nz != gridPtr_->Nz) {
        throw std::invalid_argument("Fields sizes do not match.");
    }
    for (size_t i = 0; i < data_.size(); ++i) {
        data_[i] += other.data_[i];
    }
    // todo - workaround ??: std::execution::par_unseq not available for clang 15 (libc++ 15)
    /*std::for_each(std::execution::par_unseq, data_.begin(), data_.end(),
                  [&other, this, n = size_t(0)](Scalar &elem) mutable { elem += other.data_[n++]; });*/
}

void Field::multiply(const Field::Scalar value) {
    for (auto &elem: data_) {
        elem *= value;
    }
    // todo - workaround ??: std::execution::par_unseq not available for clang 15 (libc++ 15)
    /*std::for_each(std::execution::par_unseq, data_.begin(), data_.end(),
                  [value](Scalar &elem) { elem *= value; });*/
}


// ---------------------------------------------------------------------------------------------------------------------
// VectorField class methods
// ---------------------------------------------------------------------------------------------------------------------
void VectorField::add(const Field::Scalar value) {
    component(Axis::X).add(value);
    component(Axis::Y).add(value);
    component(Axis::Z).add(value);
}

void VectorField::add(const VectorField &other) {
    component(Axis::X).add(other.component(Axis::X));
    component(Axis::Y).add(other.component(Axis::Y));
    component(Axis::Z).add(other.component(Axis::Z));
}

void VectorField::multiply(const Field::Scalar value) {
    component(Axis::X).multiply(value);
    component(Axis::Y).multiply(value);
    component(Axis::Z).multiply(value);
}