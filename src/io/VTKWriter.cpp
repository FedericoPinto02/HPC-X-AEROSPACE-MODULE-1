#include "io/VTKWriter.hpp"
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <cstdio>


VTKWriter::VTKWriter(const OutputSettings &outputSettings, const SimulationData &simData)
    : Nx_(simData.Nx), Ny_(simData.Ny), Nz_(simData.Nz),
      dx_(simData.dx), dy_(simData.dy), dz_(simData.dz),
      enabled_(outputSettings.enabled),
      outputFrequency_(outputSettings.outputFrequency),
      baseprefix_(outputSettings.dir + "/" + outputSettings.baseFilename)
{
    if (Nx_ <= 0 || Ny_ <= 0 || Nz_ <= 0)
        throw std::invalid_argument("VTKWriter: grid dimensions must be positive");
}


bool VTKWriter::write_timestep_if_needed(size_t currStep,
                                         const Field &pressure,
                                         const VectorField &velocity)
{
    // 1. Controllo Abilitazione
    if (!enabled_) {
        return false;
    }

    // 2. Controllo Frequenza
    if (currStep == 0 || (currStep % outputFrequency_ != 0)) {
        return false;
    }

    // Preparazione dei dati e del nome del file
    int step = static_cast<int>(currStep);
    char buf[256];
    std::snprintf(buf, sizeof(buf), "%s_%04d.vtk", baseprefix_.c_str(), step);
    std::string filename = std::string(buf);

    // Creazione delle copie per la scrittura (come faceva la logica esterna)
    // Questo è cruciale per prevenire data race se si usasse un thread per la scrittura
    auto pressure_ptr = std::make_shared<Field>(pressure);
    auto velocity_ptr = std::make_shared<VectorField>(velocity);

    // Scrittura
    write_legacy(filename, pressure_ptr, velocity_ptr);

    return true; // File scritto con successo
}

// Metodo write_legacy (rimane invariato)
void VTKWriter::write_legacy(const std::string &filename,
                             const std::shared_ptr<Field> &pressure,
                             const std::shared_ptr<VectorField> &velocity) const
{

    const size_t N = static_cast<size_t>(Nx_) * Ny_ * Nz_;

    std::ofstream out(filename);
    if (!out.is_open())
        throw std::runtime_error("VTKWriter::write_legacy: cannot open file " + filename);

    out << "# vtk DataFile Version 3.0\n";
    out << "NAVIER STOKES BRINKMAN SIMULATION" << "\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << Nx_ << " " << Ny_ << " " << Nz_ << "\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << dx_ << " " << dy_ << " " << dz_ << "\n";
    out << "POINT_DATA " << N << "\n";

    // SCALAR pressure
    out << "SCALARS pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    out << std::setprecision(12);

    for (int k = 0; k < Nz_; ++k) {
        for (int j = 0; j < Ny_; ++j) {
            for (int i = 0; i < Nx_; ++i) {
                out << (*pressure)(i, j, k) << "\n";
            }
        }
    }

    // VECTORS velocity
    out << "VECTORS velocity double\n";
    const auto &vx = velocity->x();
    const auto &vy = velocity->y();
    const auto &vz = velocity->z();

    for (int k = 0; k < Nz_; ++k) {
        for (int j = 0; j < Ny_; ++j) {
            for (int i = 0; i < Nx_; ++i) {
                out << vx(i,j,k) << " " << vy(i,j,k) << " " << vz(i,j,k) << "\n";
            }
        }
    }

    out.close();
}
