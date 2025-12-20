#include "numerics/SchurSolver.hpp"

void SchurSolver::preprocess() {
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
    s00 = b_[0] + c_[0] * e_left[0];
    s01 = 0.0 + c_[0] * e_right[0];
    // - Row 1 corresponds to equation for u[N-1]
    s10 = 0.0 + a_[N - 1] * e_left[n_inner - 1];
    s11 = b_[N - 1] + a_[N - 1] * e_right[n_inner - 1];
}


void SchurSolver::solveInnerSystem(const std::vector<double> &rhs, std::vector<double> &x) {
    // We need temporary slices of a, b, c corresponding to indices 1..N-2
    std::vector<double> a_in(n_inner);
    std::vector<double> b_in(n_inner);
    std::vector<double> c_in(n_inner);

    // Copy coefficients (Offset by 1 because Inner starts at index 1)
    for (size_t i = 0; i < n_inner; ++i) {
        a_in[i] = a_[i + 1];
        b_in[i] = b_[i + 1];
        c_in[i] = c_[i + 1];
    }

    // Copy RHS into x (ThomasSolver solves in-place)
    x = rhs;

    // Run Thomas algorithm
    thomas_.solve(a_in, b_in, c_in, x);
}


std::pair<double, double> SchurSolver::condenseRHS(const std::vector<double> &f) {
    /*
     * f_s_condensed = f_s - A_si inv(A_ii) b_i
     *               | inv(A_ii) b_i = y
     *               | A_ii y = b_i
     *               = f_s - A_si y
     */

    // Extract inner part of f
    std::vector<double> f_inner(n_inner);
    for (size_t i = 0; i < n_inner; ++i) {
        f_inner[i] = f[i + 1];
    }

    // Solve A_ii y = f_i
    std::vector<double> y(n_inner);
    solveInnerSystem(f_inner, y);

    // Calculate condensed f_s:
    //  f_s_condensed = f_s - A_si y
    double f_0_cond = f[0] - (c_[0] * y[0]);
    double f_N_cond = f[N - 1] - (a_[N - 1] * y[n_inner - 1]);

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
    std::vector<double> send_S = {s00, s01, s10, s11};
    std::vector<double> send_f = {f0_cond, fN_cond};

    // 3. Receive buffers (only significant on Rank 0)
    std::vector<double> recv_S;
    std::vector<double> recv_f;

    if (lineRank_ == 0) {
        recv_S.resize(NUM_LOCAL_SCHUR_ELEMS * lineNProcs_);
        recv_f.resize(NUM_LOCAL_INTERFACES * lineNProcs_);
    }

    // 4. Gather data
    MPI_Gather(send_S.data(), NUM_LOCAL_SCHUR_ELEMS, MPI_DOUBLE,
               recv_S.data(), NUM_LOCAL_SCHUR_ELEMS, MPI_DOUBLE,
               0, lineComm_);
    MPI_Gather(send_f.data(), NUM_LOCAL_INTERFACES, MPI_DOUBLE,
               recv_f.data(), NUM_LOCAL_INTERFACES, MPI_DOUBLE,
               0, lineComm_);

    // 5. Solve global Schur system to find all shared interface values for the given line
    std::vector<double> u_s_glob(lineNProcs_ + 1);
    if (lineRank_ == 0) {
        solveGlobalInterfaceSystem(recv_S, recv_f, u_s_glob);
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


void SchurSolver::solveGlobalInterfaceSystem(const std::vector<double> &S_all,
                                             const std::vector<double> &f_all,
                                             std::vector<double> &u_s_all) {
    size_t P = lineNProcs_;
    size_t M = lineNProcs_ + 1;         // Total size of global reduced system
    std::vector<double> ga(M, 0.0);
    std::vector<double> gb(M, 0.0);
    std::vector<double> gc(M, 0.0);
    std::vector<double> rhs(M, 0.0);

    // Assembly Loop
    for (size_t p = 0; p < P; ++p) {
        // Proc p contributes to Interface p (its Left) and Interface p+1 (its Right)

        // --- Contribution to Interface p (Left) ---
        // S00 adds to diagonal of Interface p
        gb[p] += S_all[4 * p + 0];
        // S01 connects Interface p to p+1 (Upper diagonal)
        gc[p] += S_all[4 * p + 1];
        // RHS contribution
        rhs[p] += f_all[2 * p + 0];

        // --- Contribution to Interface p+1 (Right) ---
        // S11 adds to diagonal of Interface p+1
        gb[p + 1] += S_all[4 * p + 3];
        // S10 connects Interface p+1 to p (Lower diagonal)
        ga[p + 1] += S_all[4 * p + 2];
        // RHS contribution
        rhs[p + 1] += f_all[2 * p + 1];
    }

    // Solve global system
    ThomasSolver globalThomas(M);
    u_s_all = rhs; // solve in-place
    globalThomas.solve(ga, gb, gc, u_s_all);
}


void SchurSolver::solveInterior(const std::vector<double> &f, std::vector<double> &u) {
    /*
     * A_ii u_i = f_i - A_is u_s
     */

    // 1. Extract inner RHS
    std::vector<double> f_inner(n_inner);
    for (size_t i = 0; i < n_inner; ++i) {
        f_inner[i] = f[i + 1];
    }

    // 2. Modify RHS with known interface values
    // - First inner node (index 0) influenced by u[0]
    f_inner[0] -= a_[1] * u[0];
    // - Last inner node (index n_inner-1) influenced by u[N-1]
    f_inner[n_inner - 1] -= c_[N - 2] * u[N - 1];

    // 3. Solve inner system for final values
    std::vector<double> u_inner(n_inner);
    solveInnerSystem(f_inner, u_inner);

    // 4. Copy back to main vector
    for (size_t i = 0; i < n_inner; ++i) {
        u[i + 1] = u_inner[i];
    }
}