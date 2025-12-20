#include "io/VTKWriter.hpp"
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <cstdio>

VTKWriter::VTKWriter(const MpiEnv &mpi, const OutputSettings &outputSettings, const SimulationData &simData)
    : mpi_(mpi), Nx_(simData.grid->Nx), Ny_(simData.grid->Ny), Nz_(simData.grid->Nz),
      dx_(simData.grid->dx), dy_(simData.grid->dy), dz_(simData.grid->dz),
      enabled_(outputSettings.enabled),
      outputFrequency_(outputSettings.outputFrequency),
      basePrefix_(outputSettings.dir + "/" + outputSettings.baseFilename)
{
    if (Nx_ <= 0 || Ny_ <= 0 || Nz_ <= 0)
        throw std::invalid_argument("VTKWriter: grid dimensions must be positive");
}

bool VTKWriter::write_timestep_if_needed(size_t currStep,
                                         const Field &pressure,
                                         const VectorField &velocity)
{
    // 1. Check if enabled
    if (!enabled_)
    {
        return false;
    }

    // 2. Check frequency
    if (currStep != 0 && (currStep % outputFrequency_ != 0))
    {
        return false;
    }

    // Prepare data and filename
    int step = static_cast<int>(currStep);
    char buf[256];
    std::snprintf(buf, sizeof(buf), "%s_%04d_%04d.vtk", basePrefix_.c_str(), step, mpi_.rank());
    std::string filename = std::string(buf);

    // Create copies for writing (as done by external logic previously)
    // This is crucial to prevent data races if a thread were used for writing
    auto pressure_ptr = std::make_shared<Field>(pressure);
    auto velocity_ptr = std::make_shared<VectorField>(velocity);

    // Perform writing
    write_legacy(filename, pressure_ptr, velocity_ptr);

    return true; // File written successfully
}

// write_legacy method (remains unchanged)
void VTKWriter::write_legacy(const std::string &filename,
                             const std::shared_ptr<Field> &pressure,
                             const std::shared_ptr<VectorField> &velocity) const
{

    const size_t N = static_cast<size_t>(Nx_) * Ny_ * Nz_;
    const auto &grid = pressure->getGrid();
    double originX = grid.i_start * grid.dx;
    double originY = grid.j_start * grid.dy;
    double originZ = grid.k_start * grid.dz;

    std::ofstream out(filename);
    if (!out.is_open())
        throw std::runtime_error("VTKWriter::write_legacy: cannot open file " + filename);

    out << "# vtk DataFile Version 3.0\n";
    out << "NAVIER STOKES BRINKMAN SIMULATION" << "\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << Nx_ << " " << Ny_ << " " << Nz_ << "\n";
    out << "ORIGIN " << originX << " " << originY << " " << originZ << "\n";
    out << "SPACING " << dx_ << " " << dy_ << " " << dz_ << "\n";
    out << "POINT_DATA " << N << "\n";

    // SCALAR pressure
    out << "SCALARS pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    out << std::setprecision(12);

    for (int k = 0; k < Nz_; ++k)
    {
        for (int j = 0; j < Ny_; ++j)
        {
            for (int i = 0; i < Nx_; ++i)
            {
                out << (*pressure)(i, j, k) << "\n";
            }
        }
    }

    // VECTORS velocity
    out << "VECTORS velocity double\n";
    const auto &vx = velocity->component(Axis::X);
    const auto &vy = velocity->component(Axis::Y);
    const auto &vz = velocity->component(Axis::Z);

    for (int k = 0; k < Nz_; ++k)
    {
        for (int j = 0; j < Ny_; ++j)
        {
            for (int i = 0; i < Nx_; ++i)
            {
                out << vx(i, j, k) << " " << vy(i, j, k) << " " << vz(i, j, k) << "\n";
            }
        }
    }

    out.close();
}