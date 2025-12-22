#include "numerics/SchurSolver.hpp"

SchurSolver::SchurSolver(const MpiEnv &env,
                         const Axis axis,
                         const TridiagMat &matrix)
        : N(matrix.getSize()), n_inner(N - 2),
          a_(matrix.getDiag(-1)), b_(matrix.getDiag(0)), c_(matrix.getDiag(1)),
          thomas_(matrix.getSize()) {
    lineComm_ = env.lineComm(axis);
    lineNProcs_ = env.lineSize(axis);
    lineRank_ = env.lineRank(axis);

    a_in_ = std::vector<double>(a_.begin() + 1, a_.end() - 1);
    b_in_ = std::vector<double>(b_.begin() + 1, b_.end() - 1);
    c_in_ = std::vector<double>(c_.begin() + 1, c_.end() - 1);
    f_in_.resize(n_inner);
    y_.resize(n_inner);

    if (N < 3) { throw std::runtime_error("Local grid too small for Schur!"); }

    if (lineRank_ == 0) {
        size_t M = lineNProcs_ + 1;         // Total size of global reduced (interface) system
        Sa_glob_.resize(M, 0.0);
        Sb_glob_.resize(M, 0.0);
        Sc_glob_.resize(M, 0.0);
    }
}


void SchurSolver::preprocess() {
    //==================================================================================================================
    // --- Local Schur Complement matrix (2x2) -------------------------------------------------------------------------
    //      S_local = [ s00  s01 ]
    //                [ s10  s11 ]
    //==================================================================================================================
    /*
     * S = A_ss - A_si inv(A_ii) A_is
     *   | -inv(A_ii) A_is = E
     *   | A_ii E = -A_is
     *       | A_ii E_left = -A_is_left     (E_left is an unknown vector of dimension n_inner X 1)
     *       | A_ii E_right = -A_is_right   (E_right is an unknown vector of dimension n_inner X 1)
     *   = A_ss + A_si E
     */

    // 0. Initialize E_left and E_right, representing E = inv(A_ii) A_is
    std::vector<double> e_left(n_inner);    // how the left interface influences the interior
    std::vector<double> e_right(n_inner);   // how the right interface influences the interior

    // 1. Construct RHS for E_left (influence of left interface u[0]):
    //  A_ii E_left = -A_is_left
    // * The term connecting u[0] to u[1] is a[1]. Moved to RHS -> -a[1].
    std::vector<double> rhs_L(n_inner, 0.0);
    rhs_L[0] = -a_[1];

    // 2. Construct RHS for E_right (influence of right interface u[N-1]):
    //  A_ii E_right = -A_is_right
    // * The term connecting u[N-1] to u[N-2] is c[N-2]. Moved to RHS -> -c[N-2].
    std::vector<double> rhs_R(n_inner, 0.0);
    rhs_R[n_inner - 1] = -c_[N - 2];

    // 3. Solve inner systems: compute E
    solveInnerSystem(rhs_L, e_left);
    solveInnerSystem(rhs_R, e_right);

    // 4. Calculate Schur complement:
    //  S = A_ss + A_si E
    // - Row 0 corresponds to equation for u[0]
    double s00 = b_[0] + c_[0] * e_left[0];
    double s01 = 0.0 + c_[0] * e_right[0];
    // - Row 1 corresponds to equation for u[N-1]
    double s10 = 0.0 + a_[N - 1] * e_left[n_inner - 1];
    double s11 = b_[N - 1] + a_[N - 1] * e_right[n_inner - 1];


    //==================================================================================================================
    // --- Global (rank 0 only) Schur matrix [(nProc+1)x(nProc+1)] -----------------------------------------------------
    //==================================================================================================================
    // 1. Prepare data for Schur system global solving (Gather to Rank 0)
    std::vector<double> send_S = {s00, s01, s10, s11};
    std::vector<double> recv_S;

    if (lineRank_ == 0) { recv_S.resize(NUM_LOCAL_SCHUR_ELEMS * lineNProcs_); }

    MPI_Gather(send_S.data(), NUM_LOCAL_SCHUR_ELEMS, MPI_DOUBLE,
               recv_S.data(), NUM_LOCAL_SCHUR_ELEMS, MPI_DOUBLE,
               0, lineComm_);

    // 2. Assemble global Schur matrix (Rank 0)
    if (lineRank_ == 0) {
        // Reset
        std::fill(Sa_glob_.begin(), Sa_glob_.end(), 0.0);
        std::fill(Sb_glob_.begin(), Sb_glob_.end(), 0.0);
        std::fill(Sc_glob_.begin(), Sc_glob_.end(), 0.0);

        // Assembly Loop
        for (size_t p = 0; p < lineNProcs_; ++p) {
            // Proc p contributes to Interface p (its Left) and Interface p+1 (its Right)
            // --- Contribution to Interface p (Left) ---
            Sb_glob_[p] += recv_S[4 * p + 0];                       // S00 adds to diagonal of Interface p
            Sc_glob_[p] += recv_S[4 * p + 1];                       // S01 connects Interface p to p+1 (Upper diagonal)
            // --- Contribution to Interface p+1 (Right) ---
            Sb_glob_[p + 1] += recv_S[4 * p + 3];                   // S11 adds to diagonal of Interface p+1
            Sa_glob_[p + 1] += recv_S[4 * p + 2];                   // S10 connects Interface p+1 to p (Lower diagonal)
        }
    }
}


void SchurSolver::solve(const std::vector<double> &f, std::vector<double> &u) {
    if (u.size() != N) { u.resize(N); }
    // Phase 1: Solve for the interfaces (global communication required)
    solveInterface(f, u);
    // Phase 2: Solve for the interior (purely local)
    solveInterior(f, u);
}


std::pair<double, double> SchurSolver::condenseRHS(const std::vector<double> &f) {
    /*
     * f_s_condensed = f_s - A_si inv(A_ii) b_i
     *               | inv(A_ii) b_i = y
     *               | A_ii y = b_i
     *               = f_s - A_si y
     */

    // Extract inner part of f
    for (size_t i = 0; i < n_inner; ++i) {
        f_in_[i] = f[i + 1];
    }

    // Solve A_ii y = f_i
    solveInnerSystem(f_in_, y_);

    // Calculate condensed f_s:
    //  f_s_condensed = f_s - A_si y
    double f_0_cond = f[0] - (c_[0] * y_[0]);
    double f_N_cond = f[N - 1] - (a_[N - 1] * y_[n_inner - 1]);

    return {f_0_cond, f_N_cond};
}


void SchurSolver::solveInterface(const std::vector<double> &f, std::vector<double> &u) {
    /*
     * S_all u_s_all = f_s_condensed_all
     * - S_all: global Schur matrix (~PxP block-tridiagonal)
     * - u_s_all: global interface unknowns (~P)
     * - f_s_condensed_all: global condensed RHS (~P)
     */

    // 1. Condense *local* RHS for the *local* contribution to the global Schur system
    auto [f0_cond, fN_cond] = condenseRHS(f);

    // 2. Prepare data for Schur system global solving (Gather to Rank 0)
    std::vector<double> send_f = {f0_cond, fN_cond};

    // 3. Receive buffers (only significant on Rank 0)
    std::vector<double> recv_f;

    if (lineRank_ == 0) {
        recv_f.resize(NUM_LOCAL_INTERFACES * lineNProcs_);
    }

    // 4. Gather data
    MPI_Gather(send_f.data(), NUM_LOCAL_INTERFACES, MPI_DOUBLE,
               recv_f.data(), NUM_LOCAL_INTERFACES, MPI_DOUBLE,
               0, lineComm_);

    // 5. Solve global Schur system to find all shared interface values for the given line
    std::vector<double> u_s_glob(lineNProcs_ + 1);
    if (lineRank_ == 0) {
        solveGlobalInterfaceSystem(recv_f, u_s_glob);
    }

    // 6. Scatter results back to local processors, to get own local interface unknowns' values
    std::vector<double> u_s_loc(NUM_LOCAL_INTERFACES);
    if (lineRank_ == 0) {
        // Send to self (Rank 0)
        u_s_loc[0] = u_s_glob[0];
        u_s_loc[1] = u_s_glob[1];

        // Send to others
        for (int p = 1; p < lineNProcs_; ++p) {
            std::vector<double> pair = {u_s_glob[p], u_s_glob[p + 1]};
            MPI_Send(pair.data(), NUM_LOCAL_INTERFACES, MPI_DOUBLE, p, 0, lineComm_);
        }
    } else {
        MPI_Recv(u_s_loc.data(), NUM_LOCAL_INTERFACES, MPI_DOUBLE, 0, 0, lineComm_, MPI_STATUS_IGNORE);
    }

    // 7. Store interface results in u (the final solution vector)
    u[0] = u_s_loc[0];
    u[N - 1] = u_s_loc[1];
}


void SchurSolver::solveGlobalInterfaceSystem(const std::vector<double> &f_all,
                                             std::vector<double> &u_s_all) {
    size_t P = lineNProcs_;
    size_t M = lineNProcs_ + 1;         // Total size of global reduced system
    std::vector<double> rhs(M, 0.0);

    // Assembly Loop
    for (size_t p = 0; p < P; ++p) {
        // Proc p contributes to Interface p (its Left) and Interface p+1 (its Right)
        rhs[p] += f_all[2 * p + 0];             // --- RHS contribution to Interface p (Left) ---
        rhs[p + 1] += f_all[2 * p + 1];         // --- RHS contribution to Interface p+1 (Right) ---
    }

    // Solve global system
    ThomasSolver globalThomas(M);
    u_s_all = rhs; // solve in-place
    globalThomas.solve(Sa_glob_, Sb_glob_, Sc_glob_, u_s_all);
}


void SchurSolver::solveInterior(const std::vector<double> &f, std::vector<double> &u) {
    /*
     * A_ii u_i = f_i - A_is u_s
     */

    // 1. Extract inner RHS
    for (size_t i = 0; i < n_inner; ++i) {
        f_in_[i] = f[i + 1];
    }

    // 2. Modify RHS with known interface values
    // - First inner node (index 0) influenced by u[0]
    f_in_[0] -= a_[1] * u[0];
    // - Last inner node (index n_inner-1) influenced by u[N-1]
    f_in_[n_inner - 1] -= c_[N - 2] * u[N - 1];

    // 3. Solve inner system for final values
    solveInnerSystem(f_in_, y_);        // y is just the scratchpad to avoid reallocation !

    // 4. Copy back to main vector
    for (size_t i = 0; i < n_inner; ++i) {
        u[i + 1] = y_[i];
    }
}