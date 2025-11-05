#include "core/Fields.hpp"
#include <stdexcept>

// ----------------------------------------------------------------
// Field class methods
// ----------------------------------------------------------------
Field::Scalar &
Field::valueWithOffset(const size_t i, const size_t j, const size_t k,
                       const Axis offsetDirection, const int offset) {
    switch (offsetDirection) {
        case Axis::X: {
            const auto i_new = static_cast<size_t>(static_cast<ssize_t>(i) + offset);
            return (*this)(i_new, j, k);
        }
        case Axis::Y: {
            const auto j_new = static_cast<size_t>(static_cast<ssize_t>(j) + offset);
            return (*this)(i, j_new, k);
        }
        case Axis::Z: {
            const auto k_new = static_cast<size_t>(static_cast<ssize_t>(k) + offset);
            return (*this)(i, j, k_new);
        }
        default:
            throw std::invalid_argument("Invalid direction.");
    }
}

const Field::Scalar &
Field::valueWithOffset(const size_t i, const size_t j, const size_t k,
                       const Axis offsetDirection, const int offset) const {
    switch (offsetDirection) {
        case Axis::X: {
            const auto i_new = static_cast<size_t>(static_cast<ssize_t>(i) + offset);
            return (*this)(i_new, j, k);
        }
        case Axis::Y: {
            const auto j_new = static_cast<size_t>(static_cast<ssize_t>(j) + offset);
            return (*this)(i, j_new, k);
        }
        case Axis::Z: {
            const auto k_new = static_cast<size_t>(static_cast<ssize_t>(k) + offset);
            return (*this)(i, j, k_new);
        }
        default:
            throw std::invalid_argument("Invalid direction.");
    }
}

void
Field::setup(std::shared_ptr<const Grid> gridPtr, std::vector<Field::Scalar> initialValues) {
    if (!gridPtr) {
        throw std::invalid_argument("Grid pointer cannot be null.");
    }
    if (initialValues.size() != gridPtr->size()) {
        throw std::invalid_argument("Initial values size does not match grid size.");
    }
    p_grid = gridPtr;
    m_v = std::move(initialValues);
}

void Field::populate(
        const std::function<double(double x, double y, double z)> &func,
        FieldOffset offset, Axis offsetAxis
) {
    auto xOffset = (offset == FieldOffset::FACE_CENTERED && offsetAxis == Axis::X)
                   ? static_cast<double>(offset) : 0.0;
    auto yOffset = (offset == FieldOffset::FACE_CENTERED && offsetAxis == Axis::Y)
                   ? static_cast<double>(offset) : 0.0;
    auto zOffset = (offset == FieldOffset::FACE_CENTERED && offsetAxis == Axis::Z)
                   ? static_cast<double>(offset) : 0.0;
    for (auto k = 0; k < p_grid->Nz; ++k) {
        for (auto j = 0; j < p_grid->Ny; ++j) {
            for (auto i = 0; i < p_grid->Nx; ++i) {
                auto x = (i + xOffset) * p_grid->dx;
                auto y = (j + yOffset) * p_grid->dy;
                auto z = (k + zOffset) * p_grid->dz;
                (*this)(i, j, k) = func(x, y, z);
            }
        }
    }
}

void Field::populate(
        double time,
        const std::function<double(double t, double x, double y, double z)> &func,
        FieldOffset offset, Axis offsetAxis
) {
    auto xOffset = (offset == FieldOffset::FACE_CENTERED && offsetAxis == Axis::X)
                   ? static_cast<double>(offset) : 0.0;
    auto yOffset = (offset == FieldOffset::FACE_CENTERED && offsetAxis == Axis::Y)
                   ? static_cast<double>(offset) : 0.0;
    auto zOffset = (offset == FieldOffset::FACE_CENTERED && offsetAxis == Axis::Z)
                   ? static_cast<double>(offset) : 0.0;
    for (auto k = 0; k < p_grid->Nz; ++k) {
        for (auto j = 0; j < p_grid->Ny; ++j) {
            for (auto i = 0; i < p_grid->Nx; ++i) {
                auto x = (i + xOffset) * p_grid->dx;
                auto y = (j + yOffset) * p_grid->dy;
                auto z = (k + zOffset) * p_grid->dz;
                (*this)(i, j, k) = func(time, x, y, z);
            }
        }
    }
}

void
Field::reset(Field::Scalar value) {
    std::fill(m_v.begin(), m_v.end(), value);
}

void
Field::update(std::vector<Field::Scalar> newV) {
    if (newV.size() != m_v.size()) {
        throw std::invalid_argument("New vector size does not match field size.");
    }
    m_v = std::move(newV);
}

void Field::add(Field::Scalar value) {
    for (auto &elem: m_v) {
        elem += value;
    }
    // todo - workaround ??: std::execution::par_unseq not available for clang 15 (libc++ 15)
    /*std::for_each(std::execution::par_unseq, m_v.begin(), m_v.end(),
                  [value](Scalar &elem) { elem += value; });*/
}

void Field::add(Field &other) {
    if (other.getGrid()->Nx != p_grid->Nx
        || other.getGrid()->Ny != p_grid->Ny
        || other.getGrid()->Nz != p_grid->Nz) {
        throw std::invalid_argument("Fields sizes do not match.");
    }
    for (size_t i = 0; i < m_v.size(); ++i) {
        m_v[i] += other.m_v[i];
    }
    // todo - workaround ??: std::execution::par_unseq not available for clang 15 (libc++ 15)
    /*std::for_each(std::execution::par_unseq, m_v.begin(), m_v.end(),
                  [&other, this, n = size_t(0)](Scalar &elem) mutable { elem += other.m_v[n++]; });*/
}

void Field::multiply(Field::Scalar value) {
    for (auto &elem: m_v) {
        elem *= value;
    }
    // todo - workaround ??: std::execution::par_unseq not available for clang 15 (libc++ 15)
    /*std::for_each(std::execution::par_unseq, m_v.begin(), m_v.end(),
                  [value](Scalar &elem) { elem *= value; });*/
}


// ----------------------------------------------------------------
// VectorField class methods
// ----------------------------------------------------------------
Field &VectorField::operator()(Axis componentDirection) {
    switch (componentDirection) {
        case Axis::X:
            return m_x;
        case Axis::Y:
            return m_y;
        case Axis::Z:
            return m_z;
        default:
            throw std::invalid_argument("Invalid direction.");
    }
}

const Field &VectorField::operator()(Axis componentDirection) const {
    switch (componentDirection) {
        case Axis::X:
            return m_x;
        case Axis::Y:
            return m_y;
        case Axis::Z:
            return m_z;
        default:
            throw std::invalid_argument("Invalid direction.");
    }
}

Field::Scalar &VectorField::operator()(const Axis componentDirection, const size_t i, const size_t j, const size_t k) {
    switch (componentDirection) {
        case Axis::X:
            return m_x(i, j, k);
        case Axis::Y:
            return m_y(i, j, k);
        case Axis::Z:
            return m_z(i, j, k);
        default:
            throw std::invalid_argument("Invalid direction.");
    }
}

const Field::Scalar &
VectorField::operator()(const Axis componentDirection, const size_t i, const size_t j, const size_t k) const {
    switch (componentDirection) {
        case Axis::X:
            return m_x(i, j, k);
        case Axis::Y:
            return m_y(i, j, k);
        case Axis::Z:
            return m_z(i, j, k);
        default:
            throw std::invalid_argument("Invalid direction.");
    }
}

void VectorField::setup(std::shared_ptr<const Grid> gridPtr,
                        std::vector<Field::Scalar> initialX,
                        std::vector<Field::Scalar> initialY,
                        std::vector<Field::Scalar> initialZ) {
    if (!gridPtr) {
        throw std::invalid_argument("Grid pointer cannot be null.");
    }
    p_grid = gridPtr;
    m_x.setup(gridPtr, std::move(initialX));
    m_y.setup(gridPtr, std::move(initialY));
    m_z.setup(gridPtr, std::move(initialZ));
}

void VectorField::update(std::vector<Field::Scalar> newX,
                         std::vector<Field::Scalar> newY,
                         std::vector<Field::Scalar> newZ) {
    m_x.update(std::move(newX));
    m_y.update(std::move(newY));
    m_z.update(std::move(newZ));
}

void VectorField::add(Field::Scalar value) {
    m_x.add(value);
    m_y.add(value);
    m_z.add(value);
}

void VectorField::add(VectorField &other) {
    m_x.add(other.x());
    m_y.add(other.y());
    m_z.add(other.z());
}

void VectorField::multiply(Field::Scalar value) {
    m_x.multiply(value);
    m_y.multiply(value);
    m_z.multiply(value);
}