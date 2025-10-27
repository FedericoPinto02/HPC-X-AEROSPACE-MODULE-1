#ifndef NSBSOLVER_FIELDS_HPP
#define NSBSOLVER_FIELDS_HPP

#include <algorithm>
#include <cstddef>
#include <execution>
#include <memory>
#include <stdexcept>
#include <vector>

#include "core/Mesh.hpp"

  enum class Axis {
        x = 0,
        y = 1,
        z = 2
    };

/**
 * @brief Class representing a scalar field defined on a 3D grid.
 */
class Field {
public:
    using Scalar = double;

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
     * @brief Read the value in the field at given position.
     * @param i x-index
     * @param j y-index
     * @param k z-index
     * @return the value in the field at position (i,j,k)
     */
    Scalar &operator()(const size_t &i, const size_t &j, const size_t &k) {
        // Assuming row-major order
        const size_t index = i + p_grid->Nx * (j + p_grid->Ny * k);
        return m_v.at(index);
    }

    /**
     * @overload
     */
    const Scalar &operator()(const size_t &i, const size_t &j,
                             const size_t &k) const {
        // Assuming row-major order
        const size_t index = i + p_grid->Nx * (j + p_grid->Ny * k);
        return m_v.at(index);
    }

    /**
   * @brief Read the value in the field at given position considering direction offset.
   * @param i x-index
   * @param j y-index
   * @param k z-index
   * @param direction direction
   * @param offset index offset
   * @return the value in the field at position (i,j,k)
   */
    Scalar &operator()(const size_t i, const size_t j, const size_t k,const Axis direction, const int offset){
        size_t i_new = i;
        size_t j_new = j;
        size_t k_new = k;

    switch (direction) {
    case Axis::x:
      i_new = static_cast<size_t>(static_cast<ssize_t>(i) + offset);
      break;
    case Axis::y:
      j_new = static_cast<size_t>(static_cast<ssize_t>(j) + offset);
      break;
    case Axis::z:
      k_new = static_cast<size_t>(static_cast<ssize_t>(k) + offset);
      break;
    default:
      throw std::invalid_argument("Invalid direction for operator()");
    }

    return (*this)(i_new, j_new, k_new);
    };

    /**
   * @overload
   */
    const Scalar &operator()(const size_t i, const size_t j, const size_t k,const Axis direction, const int offset) const{

        size_t i_new = i;
        size_t j_new = j;
        size_t k_new = k;

    switch (direction) {
    case Axis::x:
      i_new = static_cast<size_t>(static_cast<ssize_t>(i) + offset);
      break;
    case Axis::y:
      j_new = static_cast<size_t>(static_cast<ssize_t>(j) + offset);
      break;
    case Axis::z:
      k_new = static_cast<size_t>(static_cast<ssize_t>(k) + offset);
      break;
    default:
      throw std::invalid_argument("Invalid direction for operator()");
    }

    return (*this)(i_new, j_new, k_new);
    };

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
    /**
     * @brief Pointer to the grid information (dimensions and spacing).
     */
    std::shared_ptr<const Grid> p_grid;

    /**
     * @brief Vector storing the field values in a flattened, row-major indexed 1D
     * array.
     */
    std::vector<Scalar> m_v;
};

/**
 * @brief Class representing a 3D vector field defined on a 3D grid.
 */
class VectorField {
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
     * @brief Access the x-component of the vector field.
     * @return the x-component of the vector field
     */
    Field &x() { return m_x; }

    /**
     * @brief Access the x-component value of the vector field at the given position.
     * @param i x-index
     * @param j y-index
     * @param k z-index
     * @return the x-component value of the vector field at position (i,j,k)
     */
    Field::Scalar &x(const size_t &i, const size_t &j, const size_t &k) { return m_x(i, j, k); }

    /**
     * @brief Access the y-component of the vector field.
     * @return the y-component of the vector field
     */
    Field &y() { return m_y; }

    /**
     * @brief Access the y-component value of the vector field at the given position.
     * @param i x-index
     * @param j y-index
     * @param k z-index
     * @return the y-component value of the vector field at position (i,j,k)
     */
    Field::Scalar &y(const size_t &i, const size_t &j, const size_t &k) { return m_y(i, j, k); }

    /**
     * @brief Access the z-component of the vector field.
     * @return the z-component of the vector field
     */
    Field &z() { return m_z; }

    /**
     * @brief Access the z-component value of the vector field at the given position.
     * @param i x-index
     * @param j y-index
     * @param k z-index
     * @return the z-component value of the vector field at position (i,j,k)
     */
    Field::Scalar &z(const size_t &i, const size_t &j, const size_t &k) { return m_y(i, j, k); }

    /**
     * @brief Setup the vector field with initial values for each component
     * matching a given grid.
     * @param gridPtr the pointer to the grid information
     * @param initialX the initial values to populate the x-component with
     * @param initialY the initial values to populate the y-component with
     * @param initialZ the initial values to populate the z-component with
     */
    void setup(std::shared_ptr<const Grid> gridPtr, std::vector<Field::Scalar> initialX,
               std::vector<Field::Scalar> initialY,
               std::vector<Field::Scalar> initialZ);

    /**
     * @brief Update the vector field with new vectors for each component.
     * @param newX new values to update the vector field with in the x-direction
     * @param newY new values to update the vector field with in the y-direction
     * @param newZ new values to update the vector field with in the z-direction
     */
    void update(std::vector<Field::Scalar> newX, std::vector<Field::Scalar> newY,
                std::vector<Field::Scalar> newZ);

    /**
     * @brief Add a scalar value to all components of the vector field.
     * @param value the scalar value to add
     */
    void add(Field::Scalar value);

    /**
     * @brief Add another vector field to this vector field element-wise.
     * @param other the other vector field to add
     */
    void add(VectorField &other);

    /**
     * @brief Multiply all components of the vector field by a scalar value.
     * @param value the scalar value to multiply by
     */
    void multiply(Field::Scalar value);

private:
    std::shared_ptr<const Grid> p_grid;
    Field m_x;
    Field m_y;
    Field m_z;
};

#endif // NSBSOLVER_FIELDS_HPP
