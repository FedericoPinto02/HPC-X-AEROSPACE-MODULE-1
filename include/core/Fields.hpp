#ifndef FIELDS_HPP
#define FIELDS_HPP

#include <algorithm>
#include <cstddef>
#include <execution>
#include <memory>
#include <vector>

#include "core/Mesh.hpp"


/**
 * @brief Class representing a scalar field defined on a 3D grid.
 * @tparam Scalar a floating point type
 */
template<typename Scalar = double>
class Field {
    // todo - use <concepts> if we can bump to C++20 (CMake)
    static_assert(std::is_floating_point<Scalar>::value, "Field values must be floating point types.");

public:
    Field(std::shared_ptr<Grid> gridPtr, std::vector<Scalar> initialValues)
            : p_grid(std::move(gridPtr)), m_v(std::move(initialValues)) {
        if (m_v.size() != p_grid->Nx * p_grid->Ny * p_grid->Nz) {
            throw std::invalid_argument("Initial values size does not match field dimensions.");
        }
    }

    /**
     * @brief Getter for the pointer to the grid information.
     * @return the pointer to the grid information
     */
    std::shared_ptr<const Grid> getGrid() { return p_grid; }

    /**
     * @overload
     */
    std::shared_ptr<const Grid> getGrid() const { return p_grid; }

    /**
    * @brief Read the value in the field at given position.
    * @param i x-index
    * @param j y-index
    * @param k z-index
    * @return the value in the field at position (i,j,k)
    */
    Scalar &operator()(const size_t i, const size_t j, const size_t k) {
        // Assuming row-major order
        const size_t index = i + p_grid->Nx * (j + p_grid->Ny * k);
        return m_v.at(index);
    }

    /**
     * @overload
     */
    const Scalar &operator()(const size_t i, const size_t j, const size_t k) const {
        // Assuming row-major order
        const size_t index = i + p_grid->Nx * (j + p_grid->Ny * k);
        return m_v.at(index);
    }

    /**
     * @brief Reset the field values to a specified value (default is zero).
     * @param value the value to reset the field to (default is zero)
     */
    void reset(Scalar value = Scalar(0)) {
        std::fill(m_v.begin(), m_v.end(), value);
    }

    /**
     * @brief Update the field with a new vector of values.
     * @param newV new values to update the field with
     */
    void update(std::vector<Scalar> newV) {
        if (newV.size() != m_v.size()) {
            throw std::invalid_argument("New vector size does not match field size.");
        }
        m_v = std::move(newV);
    }

    /**
     * @brief Add a scalar value to all elements in the field.
     * @param value the scalar value to add
     */
    void add(Scalar value) {
        for (auto &elem: m_v) {
            elem += value;
        }
        // todo - workaround ??: std::execution::par_unseq not available for clang 15 (libc++ 15)
        /*std::for_each(std::execution::par_unseq, m_v.begin(), m_v.end(),
                      [value](Scalar &elem) { elem += value; });*/
    }

    /**
     * @brief Add another field to this field element-wise.
     * @param other the other field to add
     */
    void add(Field<Scalar> &other) {
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

    /**
     * @brief Multiply all elements in the field by a scalar value.
     * @param value the scalar value to multiply by
     */
    void multiply(Scalar value) {
        for (auto &elem: m_v) {
            elem *= value;
        }
        // todo - workaround ??: std::execution::par_unseq not available for clang 15 (libc++ 15)
        /*std::for_each(std::execution::par_unseq, m_v.begin(), m_v.end(),
                      [value](Scalar &elem) { elem *= value; });*/
    }

private:
    /**
     * @brief Pointer to the grid information (dimensions and spacing).
     */
    const std::shared_ptr<const Grid> p_grid;

    /**
     * @brief Vector storing the field values in a flattened, row-major indexed 1D array.
     */
    std::vector<Scalar> m_v;
};


/**
 * @brief Class representing a 3D vector field defined on a 3D grid.
 * @tparam Scalar a floating point type
 */
template<typename Scalar = double>
class VectorField {
    // todo - use <concepts> if we can bump to C++20 (CMake)
    static_assert(std::is_floating_point<Scalar>::value, "Field values must be floating point types.");
public:
    VectorField(std::shared_ptr<Grid> gridPtr,
                std::vector<Scalar> initialX, std::vector<Scalar> initialY, std::vector<Scalar> initialZ
    ) : m_x(gridPtr, std::move(initialX)), m_y(gridPtr, std::move(initialY)), m_z(gridPtr, std::move(initialZ)) {
        if (!p_grid) {
            throw std::invalid_argument("Grid pointer cannot be null.");
        }
    }

    /**
     * @brief Access the x-component of the vector field.
     * @return the x-component of the vector field
     */
    Field<Scalar> &x() { return m_x; }

    /**
     * @brief Access the y-component of the vector field.
     * @return the y-component of the vector field
     */
    Field<Scalar> &y() { return m_y; }

    /**
     * @brief Access the z-component of the vector field.
     * @return the z-component of the vector field
     */
    Field<Scalar> &z() { return m_z; }

    /**
     * @brief Update the vector field with new vectors for each component.
     * @param newX new values to update the vector field with in the x-direction
     * @param newY new values to update the vector field with in the y-direction
     * @param newZ new values to update the vector field with in the z-direction
     */
    void update(std::vector<Scalar> newX, std::vector<Scalar> newY, std::vector<Scalar> newZ) {
        m_x.update(std::move(newX));
        m_y.update(std::move(newY));
        m_z.update(std::move(newZ));
    }


    /**
     * @brief Add a scalar value to all components of the vector field.
     * @param value the scalar value to add
     */
    void add(Scalar value) {
        m_x.add(value);
        m_y.add(value);
        m_z.add(value);
    }

    /**
     * @brief Add another vector field to this vector field element-wise.
     * @param other the other vector field to add
     */
    void add(VectorField &other) {
        m_x.add(other.x());
        m_y.add(other.y());
        m_z.add(other.z());
    }

    /**
     * @brief Multiply all components of the vector field by a scalar value.
     * @param value the scalar value to multiply by
     */
    void multiply(Scalar value) {
        m_x.multiply(value);
        m_y.multiply(value);
        m_z.multiply(value);
    }

private:
    std::shared_ptr<const Grid> p_grid;
    Field<Scalar> m_x;
    Field<Scalar> m_y;
    Field<Scalar> m_z;
};

#endif // FIELDS_HPP
