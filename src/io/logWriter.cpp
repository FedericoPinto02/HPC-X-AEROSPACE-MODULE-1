#include "io/logWriter.hpp"
#include <iomanip>
#include <chrono>

LogWriter::LogWriter(const LoggingSettings &logSettings)
    : logToFile_(logSettings.logToFile),
      logToConsole_(logSettings.logToConsole),
      logDir_(logSettings.dir),
      filename_(logSettings.filename)
{
    if (logToFile_) {
      std::string fullPath = logDir_ + "/" + filename_;
      file_.open(fullPath, std::ios::out | std::ios::trunc);
      
      if (!file_.is_open()) {
          std::cerr << "[ERROR] Could not open log file.\n";
          logToFile_ = false;
      }
    }
}

void LogWriter::write(const std::string& msg) {
  if (logToConsole_) std::cout << msg;
  if (logToFile_ && file_.is_open()) file_ << msg << std::flush;
}

void LogWriter::printSimulationHeader(const InputData &input, const SimulationData &simData)
{
    std::ostringstream s;
    s << "\n";
    s << separator(60, '#');
    s << "   NSB SOLVER - NAVIER STOKES BRINKMAN SIMULATION   \n";
    s << separator(60, '#') << "\n";

    // SIMULATION MODE
    if (input.mesh.input_for_manufactured_solution) {
        s << "[SIMULATION MODE]\n"
          << "  Type       : MANUFACTURED SOLUTION (Validation)\n"
          << "  Note       : Grid overridden to isotropic, L = nu (Re = 1)\n\n";
    }

    // GRID INFO
    s << "[GRID CONFIGURATION]\n"
      << "  Dimensions : " << simData.Nx << " x " << simData.Ny << " x " << simData.Nz << "\n"
      << "  Domain Size: " << simData.Lx << " x " << simData.Ly << " x " << simData.Lz << "\n"
      << "  Spacing    : " << "dx=" << simData.dx << ", dy=" << simData.dy << ", dz=" << simData.dz << "\n"
      << "  Total Cells: " << (simData.Nx * simData.Ny * simData.Nz) << "\n\n";

    // TIME
    s << "[TIME SETTINGS]\n"
      << "  Time Step (dt) : " << simData.dt << "\n"
      << "  End Time       : " << simData.totalSimTime << "\n"
      << "  Total Steps    : " << simData.totalSteps << "\n\n";

    // SETUP
    s << "[SETUP]\n"
      << "  Schur Domains  : " << input.parallelization.schurDomains << "\n"
      << "  Output File    : " << input.output.baseFilename << "\n"
      << "  Log File       : " << input.logging.filename << "\n\n";
    
    s << separator() << "\n";

    write(s.str());
}

// --- TABELLA DEI PASSI ---

void LogWriter::printStepHeader() {
    std::ostringstream s;
    s << std::left 
      << std::setw(8)  << "STEP"
      << std::setw(12) << "TIME [s]" 
      << "STATUS" << "\n";
    s << separator(60, '-');
    write(s.str());
}

void LogWriter::printStepProgress(int step, double time, double dt, double elapsedSec, bool isOutputStep) 
{
    std::ostringstream s;
    
    // Usiamo fixed e setprecision per allineare i numeri
    s << std::left 
      << std::setw(8) << step
      << std::fixed << std::setprecision(3) << std::setw(12) << elapsedSec;

    if (isOutputStep) {
        s << "[VTK SAVED]";
    }
    s << "\n";
      
    write(s.str());
}