#ifndef NSBSOLVER_VTKWRITER_HPP
#define NSBSOLVER_VTKWRITER_HPP

#include <string>
#include <memory>
#include "core/Fields.hpp"

/**
 * @brief Class to write 3D fields to legacy VTK (STRUCTURED_POINTS)
 *        for a uniform Cartesian grid.
 */
class VTKWriter {
public:
    /**
     * @brief Construct a VTKWriter for a given 3D grid.
     * @param Nx number of cells in x
     * @param Ny number of cells in y
     * @param Nz number of cells in z
     * @param dx grid spacing in x
     * @param dy grid spacing in y
     * @param dz grid spacing in z
     */
    VTKWriter(int Nx, int Ny, int Nz,
              double dx, double dy, double dz);

    /**
     * @brief Write a scalar field and a vector field to a legacy VTK file.
     * @param filename path to the output file
     * @param pressure scalar field
     * @param velocity vector field
     * @param title optional title for VTK header
     */
    void write_legacy(const std::string &filename,
                      const std::shared_ptr<Field> &pressure,
                      const std::shared_ptr<VectorField> &velocity,
                      const std::string &title = "Navier-Stokes Output") const;

    /**
     * @brief Write timestep files named like baseprefix_0001.vtk
     */
    void write_timestep(const std::string &baseprefix, int step,
                        const std::shared_ptr<Field> &pressure,
                        const std::shared_ptr<VectorField> &velocity,
                        const std::string &title = "Navier-Stokes Output") const;
                        

private:
    int Nx_, Ny_, Nz_;
    double dx_, dy_, dz_;
};


#endif // NSBSOLVER_VTKWRITER_HPP
