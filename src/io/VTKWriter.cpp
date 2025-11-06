#include "io/VTKWriter.hpp"
#include "simulation/SimulationContext.hpp"
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <cstdio>


VTKWriter::VTKWriter(int Nx, int Ny, int Nz, double dx, double dy, double dz)
    : Nx_(Nx), Ny_(Ny), Nz_(Nz), dx_(dx), dy_(dy), dz_(dz)
{
    if (Nx_ <= 0 || Ny_ <= 0 || Nz_ <= 0)
        throw std::invalid_argument("VTKWriter: grid dimensions must be positive");
}

void VTKWriter::write_legacy(const std::string &filename,
                             const std::shared_ptr<Field> &pressure,
                             const std::shared_ptr<VectorField> &velocity,
                             const std::string &title) const
{

    const size_t N = static_cast<size_t>(Nx_) * Ny_ * Nz_;

    std::ofstream out(filename);
    if (!out.is_open())
        throw std::runtime_error("VTKWriter::write_legacy: cannot open file " + filename);

    out << "# vtk DataFile Version 3.0\n";
    out << title << "\n";
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

void VTKWriter::write_timestep(const std::string &baseprefix, int step,
                               const std::shared_ptr<Field> &pressure,
                               const std::shared_ptr<VectorField> &velocity,
                               const std::string &title) const
{
    char buf[256];
    std::snprintf(buf, sizeof(buf), "%s_%04d.vtk", baseprefix.c_str(), step);
    write_legacy(std::string(buf), pressure, velocity, title);
}

