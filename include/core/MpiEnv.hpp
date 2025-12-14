#pragma once

#include <array>
#include <iostream>
#include <mpi.h>


/// Enum class representing the three coordinate axes.
enum class Axis {
    X = 0, Y = 1, Z = 2
};

const size_t AXIS_COUNT = 3;

/// Class that encapsulates MPI environment and Cartesian topology management.
class MpiEnv {
private:
    /// The rank of the current process.
    int rank_;
    /// The total number of processes.
    int size_;

    /// The communicator with Cartesian topology.
    MPI_Comm cartComm_ = MPI_COMM_NULL;
    /// The dimensions of the Cartesian grid.
    std::array<int, AXIS_COUNT> dims_ = {0, 0, 0};
    /// The coordinates of the current process in the Cartesian topology.
    std::array<int, AXIS_COUNT> coords_ = {0, 0, 0};

    /// The line sub-communicators, one for each axis.
    std::array<MPI_Comm, AXIS_COUNT> lineComms_{MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL};
    /// The ranks of the current process in each line sub-communicator.
    std::array<int, AXIS_COUNT> lineRanks_{0, 0, 0};

public:
    /**
     * @brief Constructor that initializes the MPI environment.
     * @param argc the argument count from main
     * @param argv the argument vector from main
     */
    MpiEnv(int &argc, char **&argv) {
        int isInitialized;
        MPI_Initialized(&isInitialized);

        if (!isInitialized) {
            MPI_Init(&argc, &argv);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
            MPI_Comm_size(MPI_COMM_WORLD, &size_);
        }
    }

    /// Destructor that finalizes the MPI environment.
    ~MpiEnv() {
        int isFinalized;
        MPI_Finalized(&isFinalized);
        if (!isFinalized) { MPI_Finalize(); }
    }

    // Delete copy/move constructors to prevent multiple finalizations
    MpiEnv(const MpiEnv &) = delete;

    MpiEnv &operator=(const MpiEnv &) = delete;

    //==================================================================================================================
    //--- Main Communicator --------------------------------------------------------------------------------------------
    //==================================================================================================================
    /// Getter for the current process rank.
    [[nodiscard]] inline int rank() const { return rank_; }

    /// Getter for the total number of processes.
    [[nodiscard]] inline int size() const { return size_; }

    /// Getter for the Cartesian communicator.
    [[nodiscard]] inline MPI_Comm cartComm() const { return cartComm_; }

    /**
     * @brief Setup the MPI Cartesian topology with given dimensions and periodicity.
     * @param dims the dimensions of the Cartesian grid (if 0, MPI will choose the best decomposition based on the
     *  communicator size).
     * @param periods the periodicity of each dimension
     * @param reorder whether to allow MPI to reorder ranks for better locality
     */
    void setupTopology(
            std::array<int, AXIS_COUNT> dims = {0, 0, 0},
            std::array<int, AXIS_COUNT> periods = {0, 0, 0},
            int reorder = 1
    );

    /// Getter for the dimensions of the Cartesian grid.
    [[nodiscard]] inline const std::array<int, AXIS_COUNT> &dims() const { return dims_; }

    /// Getter for the coordinates of the current process in the Cartesian topology.
    [[nodiscard]] inline const std::array<int, AXIS_COUNT> &coords() const { return coords_; }

    /**
     * @brief Get the rank of the neighboring process in a given direction and displacement,
     *  in the Cartesian communicator.
     * @param direction the direction (Axis::X, Axis::Y, or Axis::Z)
     * @param disp the displacement (+1 or -1)
     * @return the rank of the neighboring process, or MPI_PROC_NULL if no neighbor exists
     */
    [[nodiscard]] int getNeighborRank(Axis direction, int disp) const;

    /// Checks if the current process is the master (rank 0).
    [[nodiscard]] inline bool isMaster() const { return rank_ == 0; }


    //==================================================================================================================
    //--- Line Communicators -------------------------------------------------------------------------------------------
    //==================================================================================================================
    /**
     * @brief Getter for the line sub-communicator along a specific axis.
     * @param axis the axis (Axis::X, Axis::Y, or Axis::Z)
     */
    [[nodiscard]] inline MPI_Comm lineComm(Axis axis) const { return lineComms_[static_cast<int>(axis)]; }

    /**
     * @brief Getter for the rank of the current process in a specific line sub-communicator.
     * @param axis the axis (Axis::X, Axis::Y, or Axis::Z)
     */
    [[nodiscard]] inline int lineRank(Axis axis) const {
        return lineRanks_[static_cast<int>(axis)];
    }

    [[nodiscard]] inline int lineSize(Axis axis) const {
        int lineSize;
        MPI_Comm_size(lineComms_[static_cast<int>(axis)], &lineSize);
        return lineSize;
    }

    /**
     * @brief Checks if the current process is the master along a specific axis (coordinate is 0).
     * @param axis the axis (Axis::X, Axis::Y, or Axis::Z)
     */
    [[nodiscard]] inline bool isLineMaster(Axis axis) const {
        return lineRanks_[static_cast<int>(axis)] == 0;
    }
};

