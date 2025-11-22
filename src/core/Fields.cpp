#include "core/Fields.hpp"

#include <stdexcept>

// ---------------------------------------------------------------------------------------------------------------------
// Field class methods
// ---------------------------------------------------------------------------------------------------------------------
Field::Scalar &
Field::valueWithOffset(size_t i, size_t j, size_t k, Axis offsetDirection, int offset) {
    switch (offsetDirection) {
        case Axis::X: {
            const auto i_new = static_cast<size_t>(static_cast<ssize_t>(i) + offset);
            return data_[idx(i_new, j, k)];
        }
        case Axis::Y: {
            const auto j_new = static_cast<size_t>(static_cast<ssize_t>(j) + offset);
            return data_[idx(i, j_new, k)];
        }
        case Axis::Z: {
            const auto k_new = static_cast<size_t>(static_cast<ssize_t>(k) + offset);
            return data_[idx(i, j, k_new)];
        }
        default:
            throw std::invalid_argument("Invalid direction.");
    }
}

const Field::Scalar &
Field::valueWithOffset(size_t i, size_t j, size_t k, Axis offsetDirection, int offset) const {
    switch (offsetDirection) {
        case Axis::X: {
            const auto i_new = static_cast<size_t>(static_cast<ssize_t>(i) + offset);
            return data_[idx(i_new, j, k)];
        }
        case Axis::Y: {
            const auto j_new = static_cast<size_t>(static_cast<ssize_t>(j) + offset);
            return data_[idx(i, j_new, k)];
        }
        case Axis::Z: {
            const auto k_new = static_cast<size_t>(static_cast<ssize_t>(k) + offset);
            return data_[idx(i, j, k_new)];
        }
        default:
            throw std::invalid_argument("Invalid direction.");
    }
}

void Field::populate(double time) {
    for (auto k = 0; k < grid_->Nz; ++k) {
        for (auto j = 0; j < grid_->Ny; ++j) {
            for (auto i = 0; i < grid_->Nx; ++i) {
                auto x = grid_->to_x(i, offset_, offsetAxis_);
                auto y = grid_->to_y(j, offset_, offsetAxis_);
                auto z = grid_->to_z(k, offset_, offsetAxis_);
                data_[idx(i, j, k)] = populateFunction_(x, y, z, time);
            }
        }
    }
}

void
Field::reset(const Field::Scalar value) {
    std::fill(data_.begin(), data_.end(), value);
}

void
Field::update(std::vector<Field::Scalar> newV) {
    if (newV.size() != data_.size()) {
        throw std::invalid_argument("New vector size does not match field size.");
    }
    data_ = std::move(newV);
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
    if (other.getGrid().Nx != grid_->Nx
        || other.getGrid().Ny != grid_->Ny
        || other.getGrid().Nz != grid_->Nz) {
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
Field::Scalar &VectorField::operator()(Axis componentDirection, size_t i, size_t j, size_t k) {
    return component(componentDirection).value(i, j, k);
}

const Field::Scalar &VectorField::operator()(Axis componentDirection, size_t i, size_t j, size_t k) const {
    return component(componentDirection).value(i, j, k);
}

void VectorField::update(std::vector<Field::Scalar> newX,
                         std::vector<Field::Scalar> newY,
                         std::vector<Field::Scalar> newZ) {
    component(Axis::X).update(std::move(newX));
    component(Axis::Y).update(std::move(newY));
    component(Axis::Z).update(std::move(newZ));
}

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