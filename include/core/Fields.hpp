#ifndef NSBSOLVER_FIELDS_HPP
#define NSBSOLVER_FIELDS_HPP

#include <algorithm>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>
#include "core/Mesh.hpp"


/**
 * @brief Enum class representing the three coordinate axes.
 */
enum class Axis {
    X = 0,
    Y = 1,
    Z = 2
};
/**
 * @brief Enum class representing the three coordinate axes.
 */
enum class Axis {
    X = 0,
    Y = 1,
    Z = 2
};

/**
 * @brief Class representing a scalar field defined on a 3D grid.
 */
class Field {
public:
    using Scalar = double;
    
    
    std::vector<Field::Scalar>& getData() { return m_v; }
    const std::vector<Field::Scalar>& getData() const { return m_v; }

    /**
     * @brief Getter for the pointer to the grid information.
     * @return the pointer to the grid information
     */
    [[nodiscard]] std::shared_ptr<const Grid> getGrid() { return p_grid; }

    /**
     * @overload
     */
    [[nodiscard]] std::shared_ptr<const Grid> getGrid() const { return p_grid; }

    /**
     * @brief Access the value in the field at given position.
     * @param i x-index
     * @param j y-index
     * @param k z-index
     * @return the value in the field at position (i,j,k)
     */
    [[nodiscard]] Scalar &operator()(const size_t i, const size_t j, const size_t k) {
        // Assuming row-major order
        const size_t index = i + p_grid->Nx * (j + p_grid->Ny * k);
        return m_v[index];
    }

    /**
     * @overload
     */
    [[nodiscard]] const Scalar &operator()(const size_t i, const size_t j, const size_t k) const {
        // Assuming row-major order
        const size_t index = i + p_grid->Nx * (j + p_grid->Ny * k);
        return m_v[index];
    }

    /**
     * @copydoc Field::operator()(size_t, size_t, size_t)
     *
     * This method is a named alias for the access operator.
     */
    [[nodiscard]] Scalar &value(const size_t i, const size_t j, const size_t k) {
        return this->operator()(i, j, k);
    }

    /**
     * @overload
     */
    [[nodiscard]] const Scalar &value(const size_t i, const size_t j, const size_t k) const {
        return this->operator()(i, j, k);
    }

    /**
     * @brief Access the value in the field at given position considering direction offset.
     * @param i the x-index
     * @param j the y-index
     * @param k the z-index
     * @param offsetDirection the direction for the offset (Axis::X, Axis::Y, or Axis::Z)
     * @param offset the offset value (can be positive or negative)
     * @return the value in the field at position (i,j,k)
     */
    [[nodiscard]] Scalar &
    valueWithOffset(const size_t i, const size_t j, const size_t k, const Axis offsetDirection, const int offset);

    /**
     * @overload
     */
    [[nodiscard]] const Scalar &
    valueWithOffset(const size_t i, const size_t j, const size_t k, const Axis offsetDirection, const int offset) const;

    /**
     * @brief Setup the field with initial values matching a given grid.
     * @param gridPtr the pointer to the grid information
     * @param initialValues the initial values to populate the field with
     */
    void setup(std::shared_ptr<const Grid> gridPtr,
               std::vector<Scalar> initialValues);

    /**
     * @brief Reset the field values to a specified value (default is zero).
     * @param value the value to reset the field to (default is zero)
     */
    void reset(Scalar value = Scalar(0));

    /**
     * @brief Update the field with a new vector of values.
     * @param newV new values to update the field with
     */
    void update(std::vector<Scalar> newV);

    /**
     * @brief Add a scalar value to all elements in the field.
     * @param value the scalar value to add
     */
    void add(Scalar value);

    /**
     * @brief Add another field to this field element-wise.
     * @param other the other field to add
     */
    void add(Field &other);

    /**
     * @brief Multiply all elements in the field by a scalar value.
     * @param value the scalar value to multiply by
     */
    void multiply(Scalar value);

private:
    /// Pointer to the grid information (dimensions and spacing).
    std::shared_ptr<const Grid> p_grid;

    /// Vector storing the field values in a flattened, row-major indexed 1D array.
    std::vector<Scalar> m_v;
};


/**
 * @brief Class representing a 3D vector field defined on a 3D grid.
 */
class VectorField {
    using Scalar = Field::Scalar;
public:
    /**
     * @brief Getter for the pointer to the grid information.
     * @return the pointer to the grid information
     */
    [[nodiscard]] std::shared_ptr<const Grid> getGrid() { return p_grid; }

    /**
     * @overload
     */
    [[nodiscard]] std::shared_ptr<const Grid> getGrid() const { return p_grid; }

    /**
     * @brief Access the component of the vector field in the specified direction.
     * @param componentDirection the direction (Axis::X, Axis::Y, or Axis::Z)
     * @return the component of the vector field in the specified direction
     */
    Field &operator()(const Axis componentDirection);

    /**
     * @overload
     */
    const Field &operator()(const Axis componentDirection) const;

    /**
     * @copydoc VectorField::operator()(Axis)
     *
     * This method is a named alias for the access operator.
     */
    [[nodiscard]] Field &component(const Axis componentDirection) {
        return this->operator()(componentDirection);
    }

    /**
     * @overload
     */
    [[nodiscard]] const Field &
    component(const Axis componentDirection) const {
        return this->operator()(componentDirection);
    }

    /**
     * @brief Access the component value of the vector field in the specified direction at the given position.
     * @param direction the direction (Axis::X, Axis::Y, or Axis::Z)
     * @param i x-index
     * @param j y-index
     * @param k z-index
     * @return the component value of the vector field in the specified direction at position (i,j,k)
     */
    Scalar &operator()(const Axis componentDirection, const size_t i, const size_t j, const size_t k);

    /**
     * @overload
     */
    const Scalar &operator()(const Axis componentDirection, const size_t i, const size_t j, const size_t k) const;

    /**
     * @copydoc VectorField::operator()(Axis, size_t, size_t, size_t)
     *
     * This method is a named alias for the access operator.
     */
    [[nodiscard]] Scalar &value(const Axis componentDirection, const size_t i, const size_t j, const size_t k) {
        return this->operator()(componentDirection, i, j, k);
    }

    /**
     * @overload
     */
    [[nodiscard]] const Scalar &
    value(const Axis componentDirection, const size_t i, const size_t j, const size_t k) const {
        return this->operator()(componentDirection, i, j, k);
    }

    /**
     * @deprecated
     * @brief Access the x-component of the vector field.
     * @return the x-component of the vector field
     */
    Field &x() { return m_x; }

    /**
     * @deprecated
     * @brief Access the x-component value of the vector field at the given position.
     * @param i x-index
     * @param j y-index
     * @param k z-index
     * @return the x-component value of the vector field at position (i,j,k)
     */
    Scalar &x(const size_t i, const size_t j, const size_t k) { return m_x(i, j, k); }

    /**
     * @deprecated
     * @brief Access the y-component of the vector field.
     * @return the y-component of the vector field
     */
    Field &y() { return m_y; }


    /**
     * @deprecated
     * @brief Access the y-component value of the vector field at the given position.
     * @param i x-index
     * @param j y-index
     * @param k z-index
     * @return the y-component value of the vector field at position (i,j,k)
     */
    Scalar &y(const size_t &i, const size_t &j, const size_t &k) { return m_y(i, j, k); }
    /**
     * @deprecated
     * @brief Access the z-component of the vector field.
     * @return the z-component of the vector field
     */
    Field &z() { return m_z; }

    /**
     * @deprecated
     * @brief Access the z-component value of the vector field at the given position.
     * @param i x-index
     * @param j y-index
     * @param k z-index
     * @return the z-component value of the vector field at position (i,j,k)
     */
    Scalar &z(const size_t i, const size_t j, const size_t k) { return m_z(i, j, k); }

    /**
     * @brief Setup the vector field with initial values for each component
     * matching a given grid.
     * @param gridPtr the pointer to the grid information
     * @param initialX the initial values to populate the x-component with
     * @param initialY the initial values to populate the y-component with
     * @param initialZ the initial values to populate the z-component with
     */
    void setup(std::shared_ptr<const Grid> gridPtr, std::vector<Scalar> initialX,
               std::vector<Scalar> initialY,
               std::vector<Scalar> initialZ);

    /**
     * @brief Update the vector field with new vectors for each component.
     * @param newX new values to update the vector field with in the x-direction
     * @param newY new values to update the vector field with in the y-direction
     * @param newZ new values to update the vector field with in the z-direction
     */
    void update(std::vector<Scalar> newX, std::vector<Scalar> newY,
                std::vector<Scalar> newZ);

    /**
     * @brief Add a scalar value to all components of the vector field.
     * @param value the scalar value to add
     */
    void add(Scalar value);

    /**
     * @brief Add another vector field to this vector field element-wise.
     * @param other the other vector field to add
     */
    void add(VectorField &other);

    /**
     * @brief Multiply all components of the vector field by a scalar value.
     * @param value the scalar value to multiply by
     */
    void multiply(Scalar value);

private:
    /// Pointer to the grid information (dimensions and spacing).
    std::shared_ptr<const Grid> p_grid;
    /// X-component of the vector field.
    Field m_x;
    /// Y-component of the vector field.
    Field m_y;
    /// Z-component of the vector field.
    Field m_z;
};

#endif // NSBSOLVER_FIELDS_HPP