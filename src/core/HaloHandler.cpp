#include "core/HaloHandler.hpp"

#include <cstring> // for std::memcpy

void HaloHandler::exchange(VectorField &vectorField) {
    exchange(vectorField(Axis::X));
    exchange(vectorField(Axis::Y));
    exchange(vectorField(Axis::Z));
}

void HaloHandler::exchange(Field &field) {
    exchangeX(field);
    exchangeY(field);
    exchangeZ(field);
}

void HaloHandler::exchangeX(Field &field) {
    const auto &grid = field.getGrid();
    const int dim_x = 0;
    if (mpi.dims()[dim_x] <= 1) return;

    MPI_Comm comm = mpi.cartComm();
    int leftProc = MPI_PROC_NULL; // X-
    int rightProc = MPI_PROC_NULL; // X+
    MPI_Cart_shift(comm, dim_x, 1, &leftProc, &rightProc);

    // Dimensions
    const size_t H = grid.n_halo;
    const size_t Nx = grid.Nx;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz_tot = grid.Nz + 2 * H;
    // X-stride is 1, so we must iterate point-by-point or gather strided vectors.
    // Buffer Size: H points * Ny_tot lines * Nz_tot planes
    const size_t bufferSize = H * Ny_tot * Nz_tot;

    sendBufLeft_.resize(bufferSize);
    sendBufRight_.resize(bufferSize);
    recvBufLeft_.resize(bufferSize);
    recvBufRight_.resize(bufferSize);

    // --- 1. PACKING ---
    // Access data raw for speed
    double *data = field.getData().data();
    const size_t strideY = grid.Nx + 2 * H; // Nx_tot
    const size_t strideZ = strideY * Ny_tot;

    // PACK SEND TO LEFT (Source: My logical index 1 to H) -> Neighbor's Right Halo
    // We skip index 0 because it is Shared.
    if (leftProc != MPI_PROC_NULL) {
        size_t ptr = 0;
        for (size_t k = 0; k < Nz_tot; ++k) {
            for (size_t j = 0; j < Ny_tot; ++j) {
                size_t rowStart = k * strideZ + j * strideY;
                for (size_t h = 0; h < H; ++h) {
                    // Src Index: HaloOffset(H) + 1 + h
                    sendBufLeft_[ptr++] = data[rowStart + H + 1 + h];
                }
            }
        }
    }

    // PACK SEND TO RIGHT (Source: My logical index Nx-H-1 to Nx-2) -> Neighbor's Left Halo
    // We skip index Nx-1 because it is Shared.
    // We want the H points immediately preceding the shared boundary.
    // End Index = Nx + H - 2. Start Index = (Nx + H - 2) - H + 1 = Nx - 1.
    // Wait, let's trace: Shared is at local logical Nx-1.
    // We want [Nx-1-H, Nx-2]. In physical memory (shifted by H): [Nx-1, Nx+H-2].
    if (rightProc != MPI_PROC_NULL) {
        size_t ptr = 0;
        for (size_t k = 0; k < Nz_tot; ++k) {
            for (size_t j = 0; j < Ny_tot; ++j) {
                size_t rowStart = k * strideZ + j * strideY;
                for (size_t h = 0; h < H; ++h) {
                    // Src Index: (H + Nx - 1 - H) + h = Nx + h - 1 ?
                    // Let's use relative logic:
                    // Shared Node is at offset (H + Nx - 1).
                    // We want H nodes before it? No, we want the *internal* nodes.
                    // Overlapping: Neighbor's u[0] is my u[N-1].
                    // Neighbor needs my u[N-2] to fill their u[-1].
                    // So I send the point at logical Nx-2.
                    // Physical index: H + Nx - 2.
                    // If H > 1, I send [H+Nx-1-H, H+Nx-2].
                    size_t idx = (H + Nx - 1 - H) + h;
                    sendBufRight_[ptr++] = data[rowStart + idx];
                }
            }
        }
    }

    // --- 2. EXCHANGE ---
    MPI_Request reqs[4];
    MPI_Irecv(recvBufLeft_.data(), (int) bufferSize, MPI_DOUBLE, leftProc, 10, comm, &reqs[0]);
    MPI_Irecv(recvBufRight_.data(), (int) bufferSize, MPI_DOUBLE, rightProc, 20, comm, &reqs[1]);
    MPI_Isend(sendBufLeft_.data(), (int) bufferSize, MPI_DOUBLE, leftProc, 20, comm, &reqs[2]);
    MPI_Isend(sendBufRight_.data(), (int) bufferSize, MPI_DOUBLE, rightProc, 10, comm, &reqs[3]);
    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

    // --- 3. UNPACKING ---
    // Unpack Left Recv -> Into my Left Halo (Indices 0..H-1)
    if (leftProc != MPI_PROC_NULL) {
        size_t ptr = 0;
        for (size_t k = 0; k < Nz_tot; ++k) {
            for (size_t j = 0; j < Ny_tot; ++j) {
                size_t rowStart = k * strideZ + j * strideY;
                for (size_t h = 0; h < H; ++h) {
                    data[rowStart + h] = recvBufLeft_[ptr++];
                }
            }
        }
    }
    // Unpack Right Recv -> Into my Right Halo (Indices H+Nx..End)
    if (rightProc != MPI_PROC_NULL) {
        size_t ptr = 0;
        for (size_t k = 0; k < Nz_tot; ++k) {
            for (size_t j = 0; j < Ny_tot; ++j) {
                size_t rowStart = k * strideZ + j * strideY;
                for (size_t h = 0; h < H; ++h) {
                    // Dest: H + Nx + h
                    // Since shared is at H+Nx-1. Next is H+Nx.
                    data[rowStart + H + Nx + h] = recvBufRight_[ptr++];
                }
            }
        }
    }
}

void HaloHandler::exchangeY(Field &field) {
    const auto &grid = field.getGrid();
    const int dim_y = 1;
    if (mpi.dims()[dim_y] <= 1) return;

    MPI_Comm comm = mpi.cartComm();
    int botProc = MPI_PROC_NULL; // Y-
    int topProc = MPI_PROC_NULL; // Y+
    MPI_Cart_shift(comm, dim_y, 1, &botProc, &topProc);

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H; // We send FULL X-width (including X-halos)
    const size_t Ny = grid.Ny;
    const size_t Nz_tot = grid.Nz + 2 * H;

    // Buffer Size: H lines * Nx_tot points * Nz_tot planes
    const size_t bufferSize = H * Nx_tot * Nz_tot;

    sendBufLeft_.resize(bufferSize);  // Buffer for Sending to Bottom (Y-)
    sendBufRight_.resize(bufferSize); // Buffer for Sending to Top (Y+)
    recvBufLeft_.resize(bufferSize);  // Buffer for Recv from Bottom
    recvBufRight_.resize(bufferSize); // Buffer for Recv from Top

    double *data = field.getData().data();
    const size_t strideY = Nx_tot;       // Distance between j and j+1
    const size_t strideZ = strideY * (grid.Ny + 2 * H);

    // --- 1. PACKING (Using MEMCPY for speed) ---
    // SEND TO BOTTOM (My logical Y=1..H) -> Neighbor's Top Halo
    // Skip j=0 (Shared).
    if (botProc != MPI_PROC_NULL) {
        size_t ptr = 0;
        for (size_t k = 0; k < Nz_tot; ++k) {
            for (size_t h = 0; h < H; ++h) {
                // Src Line: (k, H + 1 + h)
                size_t srcOffset = k * strideZ + (H + 1 + h) * strideY;
                std::memcpy(&sendBufLeft_[ptr], &data[srcOffset], Nx_tot * sizeof(double));
                ptr += Nx_tot;
            }
        }
    }
    // SEND TO TOP (My logical Y=Ny-1-H .. Ny-2) -> Neighbor's Bottom Halo
    // Skip j=Ny-1 (Shared).
    if (topProc != MPI_PROC_NULL) {
        size_t ptr = 0;
        for (size_t k = 0; k < Nz_tot; ++k) {
            for (size_t h = 0; h < H; ++h) {
                // Src Line logic same as X: Start at (H + Ny - 1 - H).
                size_t logicalY = (Ny - 1 - H) + h; // logical index
                size_t srcOffset = k * strideZ + (H + logicalY) * strideY;
                std::memcpy(&sendBufRight_[ptr], &data[srcOffset], Nx_tot * sizeof(double));
                ptr += Nx_tot;
            }
        }
    }

    // --- 2. EXCHANGE ---
    MPI_Request reqs[4];
    MPI_Irecv(recvBufLeft_.data(), (int) bufferSize, MPI_DOUBLE, botProc, 30, comm, &reqs[0]);
    MPI_Irecv(recvBufRight_.data(), (int) bufferSize, MPI_DOUBLE, topProc, 40, comm, &reqs[1]);
    MPI_Isend(sendBufLeft_.data(), (int) bufferSize, MPI_DOUBLE, botProc, 40, comm, &reqs[2]);
    MPI_Isend(sendBufRight_.data(), (int) bufferSize, MPI_DOUBLE, topProc, 30, comm, &reqs[3]);
    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

    // --- 3. UNPACKING ---
    // Unpack From Bottom (into Y-Halo 0..H-1)
    if (botProc != MPI_PROC_NULL) {
        size_t ptr = 0;
        for (size_t k = 0; k < Nz_tot; ++k) {
            for (size_t h = 0; h < H; ++h) {
                size_t dstOffset = k * strideZ + h * strideY; // Physical 0..H-1
                std::memcpy(&data[dstOffset], &recvBufLeft_[ptr], Nx_tot * sizeof(double));
                ptr += Nx_tot;
            }
        }
    }
    // Unpack From Top (into Y-Halo H+Ny..End)
    if (topProc != MPI_PROC_NULL) {
        size_t ptr = 0;
        for (size_t k = 0; k < Nz_tot; ++k) {
            for (size_t h = 0; h < H; ++h) {
                size_t dstOffset = k * strideZ + (H + Ny + h) * strideY;
                std::memcpy(&data[dstOffset], &recvBufRight_[ptr], Nx_tot * sizeof(double));
                ptr += Nx_tot;
            }
        }
    }
}

void HaloHandler::exchangeZ(Field &field) {
    const auto &grid = field.getGrid();
    const int dim_z = 2;

    if (mpi.dims()[dim_z] <= 1) return;

    MPI_Comm comm = mpi.cartComm();
    int backProc = MPI_PROC_NULL; // Z-
    int frontProc = MPI_PROC_NULL; // Z+
    MPI_Cart_shift(comm, dim_z, 1, &backProc, &frontProc);

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz = grid.Nz;

    // We send full XY-planes (including X and Y halos)
    const size_t planeSize = Nx_tot * Ny_tot;
    const size_t bufferSize = H * planeSize;

    sendBufLeft_.resize(bufferSize);
    sendBufRight_.resize(bufferSize);
    recvBufLeft_.resize(bufferSize);
    recvBufRight_.resize(bufferSize);

    double *data = field.getData().data();

    // --- 1. PACKING (Massive MEMCPY) ---
    // SEND TO BACK (Z-) (Source: Logical k=1..H) -> Neighbor's Front Halo
    if (backProc != MPI_PROC_NULL) {
        // Contiguous block start: (H + 1) * planeSize
        size_t srcOffset = (H + 1) * planeSize;
        std::memcpy(sendBufLeft_.data(), &data[srcOffset], bufferSize * sizeof(double));
    }

    // SEND TO FRONT (Z+) (Source: Logical k=Nz-1-H .. Nz-2)
    if (frontProc != MPI_PROC_NULL) {
        size_t logicalStartK = (Nz - 1 - H);
        size_t srcOffset = (H + logicalStartK) * planeSize;
        std::memcpy(sendBufRight_.data(), &data[srcOffset], bufferSize * sizeof(double));
    }

    // --- 2. EXCHANGE ---
    MPI_Request reqs[4];
    MPI_Irecv(recvBufLeft_.data(), (int) bufferSize, MPI_DOUBLE, backProc, 50, comm, &reqs[0]);
    MPI_Irecv(recvBufRight_.data(), (int) bufferSize, MPI_DOUBLE, frontProc, 60, comm, &reqs[1]);
    MPI_Isend(sendBufLeft_.data(), (int) bufferSize, MPI_DOUBLE, backProc, 60, comm, &reqs[2]);
    MPI_Isend(sendBufRight_.data(), (int) bufferSize, MPI_DOUBLE, frontProc, 50, comm, &reqs[3]);
    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

    // --- 3. UNPACKING ---
    // Unpack From Back (into Z-Halo 0..H-1)
    if (backProc != MPI_PROC_NULL) {
        std::memcpy(data, recvBufLeft_.data(), bufferSize * sizeof(double));
    }
    // Unpack From Front (into Z-Halo H+Nz..End)
    if (frontProc != MPI_PROC_NULL) {
        size_t dstOffset = (H + Nz) * planeSize;
        std::memcpy(&data[dstOffset], recvBufRight_.data(), bufferSize * sizeof(double));
    }
}