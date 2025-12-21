#pragma once

#include <string>
#include <memory>
#include "core/Fields.hpp"
#include "io/inputReader.hpp"
#include "simulation/SimulationContext.hpp"

/**
 * @brief Class to write 3D fields to legacy VTK (STRUCTURED_POINTS)
 * for a uniform Cartesian grid, managing output frequency internally.
 */
class VTKWriter {
public:
    /**
     * @brief Construct a VTKWriter with output settings and grid dimensions.
     * @param outputSettings Output settings, including frequency.
     * @param simData Simulation data containing grid dimensions (Nx, Ny, Nz, dx, dy, dz).
     */
    VTKWriter(const MpiEnv &mpi, const OutputSettings &outputSettings, const SimulationData &simData);

    /**
     * @brief Writes a timestep file if the current step matches the output frequency.
     * It takes copies of the fields internally to prevent modification during writing.
     * @param currStep The current simulation step.
     * @param pressure The scalar pressure field.
     * @param velocity The vector velocity field.
     * @return true if the file was written, false otherwise.
     */
    bool write_timestep_if_needed(size_t currStep,
                                  const VectorField &inv_porosity,
                                  const Field &pressure,
                                  const VectorField &velocity);

private:
    const MpiEnv &mpi_;

    int Nx_, Ny_, Nz_;
    double dx_, dy_, dz_;
    bool enabled_;
    int outputFrequency_;
    std::string basePrefix_;

    /**
     * @brief Internal method to write a scalar field and a vector field to a legacy VTK file.
     */
    void write_legacy(const std::string &filename,
                      const std::shared_ptr<VectorField> &inv_porosity,
                      const std::shared_ptr<Field> &pressure,
                      const std::shared_ptr<VectorField> &velocity) const;
};