#ifndef NSBSOLVER_VTKWRITER_HPP
#define NSBSOLVER_VTKWRITER_HPP

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
     * @param outputSettings Le impostazioni di output, inclusa la frequenza.
     * @param simData I dati di simulazione che contengono le dimensioni della griglia (Nx, Ny, Nz, dx, dy, dz).
     */
    VTKWriter(const OutputSettings &outputSettings, const SimulationData &simData);

    /**
     * @brief Writes a timestep file if the current step matches the output frequency.
     * It takes copies of the fields internally to prevent modification during writing.
     * @param currStep Il passo di simulazione corrente.
     * @param pressure Il campo scalare della pressione.
     * @param velocity Il campo vettoriale della velocità.
     * @return true if the file was written, false otherwise.
     */
    bool write_timestep_if_needed(size_t currStep,
                                  const Field &pressure,
                                  const VectorField &velocity);

private:
    int Nx_, Ny_, Nz_;
    double dx_, dy_, dz_;
    bool enabled_;
    int outputFrequency_;
    std::string baseprefix_;

    /**
     * @brief Internal method to write a scalar field and a vector field to a legacy VTK file.
     */
    void write_legacy(const std::string &filename,
                      const std::shared_ptr<Field> &pressure,
                      const std::shared_ptr<VectorField> &velocity) const;
};


#endif // NSBSOLVER_VTKWRITER_HPP