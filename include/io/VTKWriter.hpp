#pragma once

#include <string>
#include <memory>
#include "core/Fields.hpp"
#include "io/inputReader.hpp"
#include "simulation/SimulationContext.hpp"

/**
 * @class VTKWriter
 * @brief Handles the export of physical fields to VTK files for visualization.
 * * This class implements the export of scalar and vector fields (Pressure, Velocity, Porosity)
 * using the `STRUCTURED_POINTS` legacy format. In a parallel execution, each MPI rank 
 * generates its own file representing its local subdomain.
 */
class VTKWriter {
public:
    /**
     * @brief Constructs a VTKWriter and initializes grid geometry for the local rank.
     * @param mpi The MPI environment context.
     * @param outputSettings Settings for output directory and base filename.
     * @param simData Context containing local grid dimensions and spacing.
     * @throws std::invalid_argument If grid dimensions (Nx, Ny, Nz) are not positive.
     */
    VTKWriter(const MpiEnv &mpi, const OutputSettings &outputSettings, const SimulationData &simData);

    /**
     * @brief Evaluates and executes VTK export based on the current simulation step.
     * * This method manages the output frequency internally. It creates deep copies of the 
     * fields (via `std::shared_ptr`) before writing. 
     * * @param currStep The current time iteration index.
     * @param inv_porosity Vector field containing inverse porosity values.
     * @param pressure Scalar field for fluid pressure.
     * @param velocity Vector field for fluid velocity components.
     * @return @c true if the file was written, @c false if the step was skipped or export is disabled.
     */
    bool write_timestep_if_needed(size_t currStep,
                                  const VectorField &inv_porosity,
                                  const Field &pressure,
                                  const VectorField &velocity);

private:
    const MpiEnv &mpi_;      ///< Reference to the MPI environment.

    int Nx_, Ny_, Nz_;       ///< Local grid dimensions.
    double dx_, dy_, dz_;    ///< Grid spacing (assumed uniform).
    bool enabled_;           ///< Global toggle for VTK export.
    int outputFrequency_;    ///< Frequency of export (every N steps).
    std::string basePrefix_; ///< Base path and prefix for output files.

    /**
     * @brief Low-level method to write data to the filesystem in ASCII Legacy VTK format.
     * * The method writes the header, coordinates (Origin/Spacing), and POINT_DATA sections.
     * Porosity is derived by taking the reciprocal of @p inv_porosity.
     * * @param filename Full path of the file to be created.
     * @param inv_porosity Shared pointer to the inverse porosity vector field.
     * @param pressure Shared pointer to the pressure scalar field.
     * @param velocity Shared pointer to the velocity vector field.
     */
    void write_legacy(const std::string &filename,
                      const std::shared_ptr<VectorField> &inv_porosity,
                      const std::shared_ptr<Field> &pressure,
                      const std::shared_ptr<VectorField> &velocity) const;
};