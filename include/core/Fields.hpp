#ifndef NSBSOLVER_FIELDS_HPP
#define NSBSOLVER_FIELDS_HPP

#include <algorithm>
#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>

#include "core/Grid.hpp"

/// Type alias for a function of four double variables (x,y,z,t) returning a double.
using Func = std::function<double(double x, double y, double z, double t)>;
const Func ZERO_FUNC = [](double /*x*/, double /*y*/, double /*z*/, double /*t*/ = 0) {
    return 0.0;
};

/// Enum describing the type of boundary condition (Neumann or Dirichlet).
enum class BoundaryType {Normal, Tangent};


/**
 * @brief Class representing a scalar field defined on a 3D grid.
 */
class Field {
public:
    using Scalar = double;

private:
    /// The pointer to the grid information.
    GridPtr gridPtr_;
    /// The offset of the field in the staggered grid.
    GridStaggering offset_;
    /// The axis where to apply the offset in the staggered grid.
    Axis offsetAxis_;

    /// Vector storing the field values in a flattened, row-major indexed 1D array.
    std::vector<Scalar> data_;

    /// The function used for populating the field.
    Func populateFunction_;

public:
    /// @deprecated - Use [idx] instead.
    std::vector<Field::Scalar> &getData() { return data_; }

    /// @deprecated - Use [idx] instead.
    const std::vector<Field::Scalar> &getData() const { return data_; }

    friend class VectorField;


    //==================================================================================================================
    //--- Setup --------------------------------------------------------------------------------------------------------
    //==================================================================================================================
    /**
     * @brief Setup the field with a function to populate it (e.g., initial condition).
     * @param grid the pointer to the grid information
     * @param populateFunction the function to use for populating the field
     * @param offset the offset of the field in the staggered grid
     * @param offsetAxis the axis where to apply the offset in the staggered grid
     */
    void setup(
            const GridPtr &grid,
            const Func &populateFunction = ZERO_FUNC,
            GridStaggering offset = GridStaggering::CELL_CENTERED,
            Axis offsetAxis = Axis::X
    ) {
        gridPtr_ = grid;
        populateFunction_ = populateFunction;
        offset_ = offset;
        offsetAxis_ = offsetAxis;
        data_.resize(gridPtr_->sizeWithHalo(), Scalar(0));
    }

    /**
     * @brief Populates the underlying field based on the stored function in (x,y,z,t).
     * @param time the time affecting the stored function in (x,y,z,t)
     */
    void populate(double time = 0);


    //==================================================================================================================
    //--- Grid accessors -----------------------------------------------------------------------------------------------
    //==================================================================================================================
    /**
     * @brief Getter for the grid information.
     * @return the pointer to the grid information
     */
    [[nodiscard]] inline const Grid &getGrid() { return *gridPtr_; }

    /// @overload
    [[nodiscard]] inline const Grid &getGrid() const { return *gridPtr_; }

    /**
     * @brief Getter for the field offset in the staggered grid.
     * @return the field offset in the staggered grid
     */
    [[nodiscard]] inline const GridStaggering &getOffset() { return offset_; }

    /// Overload
    [[nodiscard]] inline const GridStaggering &getOffset() const { return offset_; }

    /**
     * @brief Getter for the axis where the offset is applied in the staggered grid.
     * @return the axis where the offset is applied in the staggered grid
     */
    [[nodiscard]] inline const Axis &getOffsetAxis() { return offsetAxis_; }

    /// Overload
    [[nodiscard]] inline const Axis &getOffsetAxis() const { return offsetAxis_; }

    //==================================================================================================================
    //--- Data accessors -----------------------------------------------------------------------------------------------
    //==================================================================================================================
    /**
     * @brief Computes the linear index for 3D access in row-major order.
     * @param i the x-index
     * @param j the y-index
     * @param k the z-index
     * @return the corresponding 1D index
     */
    [[nodiscard]] inline size_t idx(long i, long j, long k) const {
        const size_t H = gridPtr_->n_halo;
        const size_t Nx_tot = gridPtr_->Nx + 2 * H;
        const size_t Ny_tot = gridPtr_->Ny + 2 * H;
        return ((k + H) * Ny_tot + (j + H)) * Nx_tot + (i + H);
    }

    /**
     * @brief Access the value in the field at given linear index.
     * @param index the linear index
     * @return the value in the field at the given index
     */
    [[nodiscard]] inline Scalar &operator[](size_t index) {
        return data_[index];
    }

    /// @overload
    [[nodiscard]] inline const Scalar &operator[](size_t index) const {
        return data_[index];
    }

    /**
     * @brief Access the value in the field at given position.
     * @param i the x-index
     * @param j the y-index
     * @param k the z-index
     * @return the value in the field at position (i,j,k)
     */
    [[nodiscard]] inline Scalar &operator()(long i, long j, long k) {
        return data_[idx(i, j, k)];
    }

    /// @overload
    [[nodiscard]] inline const Scalar &operator()(long i, long j, long k) const {
        return data_[idx(i, j, k)];
    }

    /**
     * @copydoc Field::operator()(long, long, long)
     *
     * This method is a named alias for the access operator.
     */
    [[nodiscard]] inline Scalar &value(long i, long j, long k) {
        return data_[idx(i, j, k)];
    }

    /// @overload
    [[nodiscard]] inline const Scalar &value(long i, long j, long k) const {
        return data_[idx(i, j, k)];
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
    [[nodiscard]] Scalar &valueWithOffset(long i, long j, long k, Axis offsetDirection, int offset);

    /// @overload
    [[nodiscard]] const Scalar &valueWithOffset(long i, long j, long k, Axis offsetDirection, int offset) const;


    //==================================================================================================================
    //--- Operations ---------------------------------------------------------------------------------------------------
    //==================================================================================================================
    /**
     * @brief Reset the field values to a specified value (default is zero).
     * @param value the value to reset the field to (default is zero)
     */
    void reset(Scalar value = Scalar(0));

    /**
     * @brief Add a scalar value to all elements in the field.
     * @param value the scalar value to add
     */
    void add(Scalar value);

    /**
     * @brief Add another field to this field element-wise.
     * @param other the other field to add
     */
    void add(const Field &other);

    /**
     * @brief Multiply all elements in the field by a scalar value.
     * @param value the scalar value to multiply by
     */
    void multiply(Scalar value);
};


/**
 * @brief Class representing a 3D vector field defined on a 3D grid.
 */
class VectorField {
public:
    using Scalar = Field::Scalar;

private:
    /// The pointer to the grid information (dimensions and spacing).
    GridPtr gridPtr_;
    /// The components of the vector field (x, y, z).
    std::array<Field, AXIS_COUNT> components_;

public:
    //==================================================================================================================
    //--- Setup --------------------------------------------------------------------------------------------------------
    //==================================================================================================================
    /**
     * @brief Setup the vector field with functions to populate each component (e.g., initial conditions).
     * @param grid the pointer to the grid information
     * @param populateXFunction the function to use for populating the x-component of the vector field
     * @param populateYFunction the function to use for populating the y-component of the vector field
     * @param populateZFunction the function to use for populating the z-component of the vector field
     */
    void setup(const GridPtr &grid,
               const Func &populateXFunction = ZERO_FUNC,
               const Func &populateYFunction = ZERO_FUNC,
               const Func &populateZFunction = ZERO_FUNC) {
        gridPtr_ = grid;
        component(Axis::X).setup(grid, populateXFunction, GridStaggering::FACE_CENTERED, Axis::X);
        component(Axis::Y).setup(grid, populateYFunction, GridStaggering::FACE_CENTERED, Axis::Y);
        component(Axis::Z).setup(grid, populateZFunction, GridStaggering::FACE_CENTERED, Axis::Z);
    }

    /**
     * @brief Populates each component of the underlying vector field based on the stored functions in (x,y,z,t).
     * @param time the time affecting the stored functions in (x,y,z,t)
     */
    void populate(double time = 0) {
        component(Axis::X).populate(time);
        component(Axis::Y).populate(time);
        component(Axis::Z).populate(time);
    }


    //==================================================================================================================
    //--- Grid accessors -----------------------------------------------------------------------------------------------
    //==================================================================================================================
    /**
     * @brief Getter for the grid information.
     * @return the pointer to the grid information
     */
    [[nodiscard]] inline const Grid &getGrid() { return *gridPtr_; }

    /// @overload
    [[nodiscard]] inline const Grid &getGrid() const { return *gridPtr_; }


    //==================================================================================================================
    //--- Data accessors -----------------------------------------------------------------------------------------------
    //==================================================================================================================
    /**
     * @brief Access the component of the vector field in the specified direction.
     * @param componentDirection the direction (Axis::X, Axis::Y, or Axis::Z)
     * @return the component of the vector field in the specified direction
     */
    [[nodiscard]] inline Field &operator()(Axis componentDirection) {
        return components_[static_cast<size_t>(componentDirection)];
    }

    /// @overload
    [[nodiscard]] inline const Field &operator()(Axis componentDirection) const {
        return components_[static_cast<size_t>(componentDirection)];
    }

    /**
     * @copydoc VectorField::operator()(Axis)
     *
     * This method is a named alias for the access operator.
     */
    [[nodiscard]] inline Field &component(Axis componentDirection) {
        return components_[static_cast<size_t>(componentDirection)];
    }

    /// @overload
    [[nodiscard]] inline const Field &component(Axis componentDirection) const {
        return components_[static_cast<size_t>(componentDirection)];
    }

    /**
     * @brief Access the component value of the vector field in the specified direction at the given position.
     * @param direction the direction (Axis::X, Axis::Y, or Axis::Z)
     * @param i the x-index
     * @param j the y-index
     * @param k the z-index
     * @return the component value of the vector field in the specified direction at position (i,j,k)
     */
    [[nodiscard]] Scalar &operator()(Axis componentDirection, long i, long j, long k);

    /// @overload
    [[nodiscard]] const Scalar &operator()(Axis componentDirection, long i, long j, long k) const;

    /**
     * @copydoc VectorField::operator()(Axis, long, long, long)
     *
     * This method is a named alias for the access operator.
     */
    [[nodiscard]] Scalar &value(Axis componentDirection, long i, long j, long k) {
        return this->operator()(componentDirection, i, j, k);
    }

    /// @overload
    [[nodiscard]] const Scalar &value(Axis componentDirection, long i, long j, long k) const {
        return this->operator()(componentDirection, i, j, k);
    }


    //==================================================================================================================
    //--- Operations ---------------------------------------------------------------------------------------------------
    //==================================================================================================================
    /**
     * @brief Add a scalar value to all components of the vector field.
     * @param value the scalar value to add
     */
    void add(Scalar value);

    /**
     * @brief Add another vector field to this vector field element-wise.
     * @param other the other vector field to add
     */
    void add(const VectorField &other);

    /**
     * @brief Multiply all components of the vector field by a scalar value.
     * @param value the scalar value to multiply by
     */
    void multiply(Scalar value);
};

#endif // NSBSOLVER_FIELDS_HPP