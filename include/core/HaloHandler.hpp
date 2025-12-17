#pragma once

#include <mpi.h>
#include <vector>

#include "core/Fields.hpp"
#include "core/MpiEnv.hpp"

/**
 * @brief Class handling halo exchanges for <code>Field</code> and <code>VectorField</code> classes in an overlapping grid.
 *
 * In an overlapping grid, the node at the boundaries are shared.
 * Therefore, the ghost cell (halo) must be filled by the neighbor's *internal* node, not the neighbor's boundary node.
 * Ignoring <code>Field</code>'s halo padding in indices:
 * <ul>
 *  <li>
 *      Left Halo (u[-1]) comes from LeftNeighbor's u[N-2] (not N-1).
 *  </li>
 *  <li>
 *      Right Halo (u[N]) comes from RightNeighbor's u[1] (not 0).
 *  </li>
 * </ul>
 */
class HaloHandler {
private:
    const MpiEnv &mpi;

    // Temporary buffers for packing/unpacking
    std::vector<double> sendBufLeft_;
    std::vector<double> sendBufRight_;
    std::vector<double> recvBufLeft_;
    std::vector<double> recvBufRight_;

public:
    explicit HaloHandler(const MpiEnv &env) : mpi(env) {}

    /**
     * @brief Update the halos for a <code>Field</code> (X, Y, then Z).
     * This ensures corners are filled correctly via sequential exchange.
     * @param field the field to exchange halos for
     */
    void exchange(Field &field);

    /**
     * @brief Update the halos for a <code>VectorField</code> (X, Y, then Z).
     * @param vectorField the vector field to exchange halos for
     */
    void exchange(VectorField &vectorField);

    /**
     * @brief Exchange halos in the X direction for a given field.
     * @param field the field to exchange halos for
     */
    void exchangeX(Field &field);

    /**
     * @brief Exchange halos in the Y direction for a given field.
     * @param field the field to exchange halos for
     */
    void exchangeY(Field &field);

    /**
     * @brief Exchange halos in the Z direction for a given field.
     * @param field the field to exchange halos for
     */
    void exchangeZ(Field &field);
};
