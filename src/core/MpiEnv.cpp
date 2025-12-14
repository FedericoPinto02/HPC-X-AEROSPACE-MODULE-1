#include "core/MpiEnv.hpp"

void MpiEnv::setupTopology(std::array<int, AXIS_COUNT> dims, std::array<int, AXIS_COUNT> periods, int reorder) {
    dims_ = dims;
    MPI_Dims_create(size_, AXIS_COUNT, dims_.data());

    if (MPI_Cart_create(MPI_COMM_WORLD,
                        AXIS_COUNT, dims_.data(),
                        periods.data(),
                        reorder, &cartComm_) != MPI_SUCCESS) {
        throw std::runtime_error("Failed to create MPI Cartesian Topology");
    }

    // Update current process rank and coordinates based on the NEW cartesian communicator
    MPI_Comm_rank(cartComm_, &rank_);
    MPI_Cart_coords(cartComm_, rank_, AXIS_COUNT, coords_.data());

    // Create a line sub-communicator for each axis
    // Store local ranks within each line communicator
    for (int axis = 0; axis < AXIS_COUNT; ++axis) {
        std::array<int, AXIS_COUNT> remain_dims = {0, 0, 0};
        remain_dims[axis] = 1; // Keep this dimension

        if (MPI_Cart_sub(cartComm_, remain_dims.data(), &lineComms_[axis])
            != MPI_SUCCESS) {
            throw std::runtime_error("Failed to create MPI line sub-communicator");
        }

        MPI_Comm_rank(lineComms_[axis], &lineRanks_[axis]);
    }
}

int MpiEnv::getNeighborRank(Axis direction, int disp) const {
    int source, dest;
    MPI_Cart_shift(cartComm_, static_cast<int>(direction), disp, &source, &dest);
    return dest;
}