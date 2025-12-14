#include "io/logWriter.hpp"
#include <iomanip>
#include <chrono>

LogWriter::LogWriter(const LoggingSettings &logSettings)
  : logToFile_(logSettings.logToFile),
    logToConsole_(logSettings.logToConsole),
    logDir_(logSettings.dir),
    filename_(logSettings.filename)
{
  if (logToFile_)
  {
    std::string fullPath = logDir_ + "/" + filename_;
    file_.open(fullPath, std::ios::out | std::ios::trunc);

    if (!file_.is_open())
    {
      std::cerr << "[ERROR] Could not open log file.\n";
      logToFile_ = false;
    }
  }
}

void LogWriter::write(const std::string &msg)
{
  if (logToConsole_)
    std::cout << msg;
  if (logToFile_ && file_.is_open())
    file_ << msg << std::flush;
}

void LogWriter::printSimulationHeader(const InputData &input, const SimulationData &simData, bool vtkWritten)
{
  std::ostringstream s;
  s << "\n";
  s << separator(60, '#');
  s << "   NSB SOLVER - NAVIER STOKES BRINKMAN SIMULATION   \n";
  s << separator(60, '#') << "\n";

  // SIMULATION MODE
  if (input.mesh.input_for_manufactured_solution)
  {
    s << "[SIMULATION MODE]\n"
      << "  Type            : MANUFACTURED SOLUTION (Validation)\n"
      << "  Note            : Grid overridden to isotropic, L = nu (Re = 1)\n\n";
  }

  // GRID INFO
  s << "[GRID CONFIGURATION]\n"
    << "  Dimensions      : " << simData.grid->Nx << " x " << simData.grid->Ny << " x " << simData.grid->Nz << "\n"
    << "  Domain Size     : " << simData.grid->Lx_glob << " x " << simData.grid->Ly_glob << " x " << simData.grid->Lz_glob << "\n"
    << "  Spacing         : " << "dx=" << simData.grid->dx << ", dy=" << simData.grid->dy << ", dz=" << simData.grid->dz << "\n"
    << "  Total Cells     : " << (simData.grid->Nx * simData.grid->Ny * simData.grid->Nz) << "\n\n";

  // TIME
  s << "[TIME SETTINGS]\n"
    << "  Time Step (dt)  : " << simData.dt << "\n"
    << "  End Time        : " << simData.totalSimTime << "\n"
    << "  Total Steps     : " << simData.totalSteps << "\n\n";

  // SETUP
  s << "[SETUP]\n"
    << "  Schur Domains   : " << input.parallelization.schurDomains << "\n"
    << "  Log File        : " << input.logging.filename << "\n\n";

  // VTK OUTPUT
  s << "[VTK OUTPUT]\n"
    << "  Enabled         : " << (input.output.enabled ? "YES" : "NO") << "\n"
    << "  Output File     : " << input.output.baseFilename << "\n"
    << "  Output Freq.    : " << input.output.outputFrequency << " timesteps\n"
    << "  First VTK Saved : " << (vtkWritten ? "YES" : "NO") << "\n\n";

  s << separator() << "\n";

  write(s.str());
}

// --- STEP TABLE ---

void LogWriter::printStepHeader()
{
  std::ostringstream s;
  s << std::left
    << std::setw(8) << "STEP"
    << std::setw(12) << "TIME"
    << std::setw(16) << "CPU TIME [s]"
    << "STATUS" << "\n";
  s << separator(60, '-');
  write(s.str());
}

void LogWriter::printStepProgress(int step, double time, double dt, double elapsedSec, bool isOutputStep)
{
  std::ostringstream s;

  // Use fixed and setprecision to align numbers
  s << std::left
    << std::setw(8) << step
    << std::fixed << std::setprecision(4) << std::setw(12) << time
    << std::fixed << std::setprecision(4) << std::setw(16) << elapsedSec;

  if (isOutputStep)
  {
    s << "[VTK SAVED]";
  }
  s << "\n";

  write(s.str());
}

void LogWriter::printFinalSummary(
  double totalCpuTimeSec,
  double meanCpuTimePerCellTimestep,
  unsigned int totalSteps,
  const unsigned int totalCells)
{
  std::ostringstream s;

  s << "\n"
    << separator(60, '=');
  s << "   SIMULATION SUMMARY & PERFORMANCE METRICS   \n";
  s << separator(60, '=') << "\n";

  s << std::fixed << std::setprecision(6);

  // I. General data
  s << "[GENERAL STATS]\n"
    << std::setw(30) << std::left << "Total Timesteps Completed:" << totalSteps << "\n"
    << std::setw(30) << std::left << "Total Grid Cells (N_cells):" << totalCells << "\n\n";

  // II. Timing
  s << "[PERFORMANCE TIMING]\n"
    << std::setw(30) << std::left << "Total CPU Time (solve loop):" << totalCpuTimeSec << " s\n"
    << std::setw(30) << std::left << "Avg CPU Time per Timestep:" << (totalCpuTimeSec / totalSteps) << " s\n";

  // III. Key metric
  s << "\n"
    << separator(60, '*');
  s << std::scientific << std::setprecision(4); // Use scientific notation for the metric

  s << "[TARGET METRIC]\n"
    << std::setw(30) << std::left << "CPU Time / (Steps * Cells):"
    << meanCpuTimePerCellTimestep << " s\n";

  s << separator(60, '*') << "\n";

  write(s.str());
}