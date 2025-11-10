#include "io/logWriter.hpp"
#include <iomanip>
#include <chrono>

LogWriter::LogWriter(const LoggingSettings &logSettings)
    : logToFile_(logSettings.logToFile),
      logToConsole_(logSettings.logToConsole),
      filename_(logSettings.filename)
{
    if (logToFile_) {
        file_.open(filename_, std::ios::out | std::ios::trunc);
        if (!file_.is_open()) {
            std::cerr << "[ERROR] Could not open log file: " << filename_ << "\n";
            logToFile_ = false;
        }
    }
}

void LogWriter::write(const std::string& msg) {
    if (logToConsole_)
        std::cout << msg;
    if (logToFile_ && file_.is_open())
        file_ << msg << std::flush;
}

void LogWriter::printInputSummary(const InputData &input)
{

    std::ostringstream s;
    s << "\n=== INPUT CONFIG SUMMARY ===\n" << std::boolalpha;

    s << "[MESH]\n"
      << "  - Grid size (Nx, Ny, Nz): " << input.mesh.nx << ", " << input.mesh.ny << ", " << input.mesh.nz << "\n"
      << "  - Spacing (dx, dy, dz):   " << input.mesh.dx << ", " << input.mesh.dy << ", " << input.mesh.dz << "\n";

    s << "[TIME]\n"
      << "  - Time step (dt):   " << input.time.dt << "\n"
      << "  - End time (t_end): " << input.time.t_end << "\n";

    s << "[PHYSICS]\n"
      << "  - Viscosity (nu):   " << input.physics.nu << "\n"
      << "  - Permeability (k): \"" << input.physics.k_expr << "\"\n";

    s << "[FORCES]\n"
      << "  - fx: \"" << input.forces.fx_expr << "\"\n"
      << "  - fy: \"" << input.forces.fy_expr << "\"\n"
      << "  - fz: \"" << input.forces.fz_expr << "\"\n";

    s << "[INITIAL CONDITIONS]\n"
      << "  - u: \"" << input.initial_conditions.u_expr << "\"\n"
      << "  - v: \"" << input.initial_conditions.v_expr << "\"\n"
      << "  - w: \"" << input.initial_conditions.w_expr << "\"\n"
      << "  - p: \"" << input.initial_conditions.p_expr << "\"\n";

    s << "[BOUNDARY CONDITIONS]\n"
      << "  - u: \"" << input.boundary_conditions.u_expr << "\"\n"
      << "  - v: \"" << input.boundary_conditions.v_expr << "\"\n"
      << "  - w: \"" << input.boundary_conditions.w_expr << "\"\n";

    s << "[OUTPUT]\n"
      << "  - Enabled:   " << input.output.enabled << "\n"
      << "  - Directory: \"" << input.output.dir << "\"\n"
      << "  - Filename:  \"" << input.output.baseFilename << "\"\n"
      << "  - Frequency: " << input.output.outputFrequency << "\n";

    s << "[LOGGING]\n"
      << "  - To File:   " << input.logging.logToFile << "\n"
      << "  - To Console:" << input.logging.logToConsole << "\n"
      << "  - Filename:  \"" << input.logging.filename << "\"\n"
      << "  - Frequency: " << input.logging.loggingFrequency << "\n";

    s << "============================\n";

    write(s.str());
}

void LogWriter::printRuntimeSummary(const SimulationData &simData)
{
    std::ostringstream s;
    s << "=== SIMULATION RUNTIME SETTINGS ===\n"
      << "[GRID]\n"
      << "  - Size: " << simData.Nx << " x " << simData.Ny << " x " << simData.Nz << "\n"
      << "  - Domain: " << simData.Lx << ", " << simData.Ly << ", " << simData.Lz << "\n"
      << "  - Spacing: " << simData.dx << ", " << simData.dy << ", " << simData.dz << "\n"
      << "[TIME]\n"
      << "  - dt: " << simData.dt << "\n"
      << "  - Total time: " << simData.totalSimTime << "\n"
      << "  - Steps: " << simData.totalSteps << "\n"
      << "[PHYSICS]\n"
      << "  - nu: " << simData.nu << "\n";

    write(s.str());
}

void LogWriter::printLoopProgress(const SimulationData& simData, double elapsedSec) 
{
  std::ostringstream s;
  s << "Step " << simData.currStep 
    << " done in " << std::fixed << std::setprecision(3) 
    << elapsedSec << " s.\n";
      
  write(s.str());
};

LogWriter::~LogWriter() {
    if (file_.is_open())
        file_.close();
}
