#include "numerics/SchurSequentialSolver.hpp"
#include <stdexcept>
#include <numeric> 
#include <cmath>   
#include <iostream> 

// --- CONSTRUCTOR (MODIFIED) ---
SchurSequentialSolver::SchurSequentialSolver(int globalSize, int numDomains, BoundaryType type)
    : nGlobal(globalSize),
      num_domains(numDomains),
      num_interfaces(numDomains > 0 ? numDomains - 1 : 0),
      bType(type)
      // schurSolver is no longer initialized here in the initializer list
{
    if (numDomains <= 0) {
        throw std::invalid_argument("Il numero di domini deve essere > 0.");
    }
    if (nGlobal < numDomains) {
        throw std::invalid_argument("La dimensione globale non può essere inferiore al numero di domini.");
    }

    // --- Initialization logic reorganized ---

    // Case P = 1 (0 interfaces)
    if (numDomains == 1) {
        domain_sizes.push_back(nGlobal);
        domain_starts.push_back(0);
        localSolvers.emplace_back(nGlobal, bType);
        // schurSolver remains nullptr, not needed
        return; 
    }

    // --- Cases P > 1 (at least 1 interface) ---

    // Initialize schurSolver ONLY if it's a system with N >= 2
    if (num_interfaces >= 2) { // P >= 3
        schurSolver = std::make_unique<LinearSys>(num_interfaces, type);
    }
    // If P=2 (num_interfaces=1), schurSolver remains nullptr.
    // We'll use schurScalarS.

    // Resize vectors
    domain_sizes.resize(num_domains);
    domain_starts.resize(num_domains);
    interface_indices.resize(num_interfaces);
    A_ie_left.resize(num_domains, 0.0);
    A_ie_right.resize(num_domains, 0.0);
    A_ei_left.resize(num_interfaces, 0.0);
    A_ei_right.resize(num_interfaces, 0.0);

    // Decomposition logic (unchanged, but now applies only for P > 1)
    int total_internal_nodes = nGlobal - num_interfaces;
    int base_internal_size = total_internal_nodes / num_domains;
    int remainder = total_internal_nodes % num_domains;

    int current_global_idx = 0;
    for (int p = 0; p < num_domains; ++p) {
        int n_i = base_internal_size + (p < remainder ? 1 : 0);
        domain_sizes[p] = n_i;
        domain_starts[p] = current_global_idx;
        
        localSolvers.emplace_back(n_i, bType); // Create local solver
        
        current_global_idx += n_i; 

        if (p < num_interfaces) {
            interface_indices[p] = current_global_idx;
            current_global_idx++; 
        }
    }
}

void SchurSequentialSolver::PreProcess(const TridiagMat& A_global) {
    if (A_global.getSize() != nGlobal) {
        throw std::runtime_error("Dimensione A_global non corrisponde a nGlobal.");
    }

    // --- Case P = 1 ---
    if (num_domains == 1) {
        localSolvers[0].matA.fillMat(
            A_global.getDiag(0), A_global.getDiag(-1), A_global.getDiag(1)
        );
        return;
    }

    // --- Cases P > 1 ---
    
    // --- Phase 1a: Extract A_ii and coefficients A_ie ---
    for (int p = 0; p < num_domains; ++p) {
        int n_i = domain_sizes[p];
        int start = domain_starts[p];

        std::vector<double> diag_i(n_i), sub_i(n_i - 1), sup_i(n_i - 1);
        for (int i = 0; i < n_i; ++i) {
            diag_i[i] = A_global.diag[start + i];
            if (i < n_i - 1) {
                sub_i[i] = A_global.subdiag[start + i];
                sup_i[i] = A_global.supdiag[start + i];
            }
        }
        localSolvers[p].matA.fillMat(diag_i, sub_i, sup_i);
        
        if (p > 0) A_ie_left[p] = A_global.subdiag[start - 1];
        if (p < num_interfaces) A_ie_right[p] = A_global.supdiag[start + n_i - 1];
    }

    // --- Phase 1b: Extract A_ei and diagonal of A_ee ---
    std::vector<double> S_diag_ee(num_interfaces); // Diagonal of A_ee
    for (int k = 0; k < num_interfaces; ++k) {
        int idx = interface_indices[k];
        A_ei_left[k]  = A_global.subdiag[idx - 1]; // a[idx]
        A_ei_right[k] = A_global.supdiag[idx]; // c[idx]
        S_diag_ee[k] = A_global.diag[idx]; // b[idx]
    }
    
    // --- Fase 2: Costruire S = A_ee - A_ei * A_ii^-1 * A_ie ---
    
    // 2a. Calcola tutti i w_l e w_r (le solve locali)
    std::vector<std::vector<double>> all_w_l(num_domains), all_w_r(num_domains);
    for (int p = 0; p < num_domains; ++p) {
        LinearSys& solver = localSolvers[p];
        int n_i = domain_sizes[p];
        
        std::vector<double> v_sparse_l = createSparseVector(n_i, 0, A_ie_left[p]);
        solver.setRhs(v_sparse_l);
        solver.ThomaSolver();
        all_w_l[p] = solver.getSolution();

        std::vector<double> v_sparse_r = createSparseVector(n_i, n_i - 1, A_ie_right[p]);
        solver.setRhs(v_sparse_r);
        solver.ThomaSolver();
        all_w_r[p] = solver.getSolution();
    }

    // 2b. Assembla S
    if (schurSolver) { // Caso P >= 3 (n_interfaces >= 2)
        std::vector<double> S_diag = S_diag_ee; // Inizia con S = A_ee
        std::vector<double> S_sub(num_interfaces - 1, 0.0);
        std::vector<double> S_sup(num_interfaces - 1, 0.0);

        // --- QUESTA È LA LOGICA CORRETTA ---
        for (int k = 0; k < num_interfaces; ++k) {
            // k è l'indice dell'interfaccia.
            // L'interfaccia 'k' è tra il dominio 'k' (sinistra) e 'k+1' (destra).
            
            // Diagonale S(k, k) = A_ee(k) 
            //   - Contributo dal dominio k (a sinistra)
            //   - Contributo dal dominio k+1 (a destra)

            // Contributo da sinistra (dominio k): A_ei_right[k] * w_r[k].back()
            S_diag[k] -= A_ei_right[k] * all_w_r[k].back();
            
            // Contributo da destra (dominio k+1): A_ei_left[k] * w_l[k+1].front()
            S_diag[k] -= A_ei_left[k] * all_w_l[k+1].front();

            // Fuori-diagonale S(k, k+1) (Sopradiagonale)
            // Questo termine collega l'interfaccia k e k+1 tramite il dominio k+1
            if (k < num_interfaces - 1) {
                int p = k + 1; // Dominio k+1
                // S(k, k+1) = -A_ei_left[k] * w_r[p].front()
                S_sup[k] = -A_ei_left[k] * all_w_r[p].front();
            }

            // Fuori-diagonale S(k, k-1) (Sottodiagonale)
            // Questo termine collega l'interfaccia k e k-1 tramite il dominio k
            if (k > 0) {
                int p = k; // Dominio k
                // S(k, k-1) = -A_ei_right[k] * w_l[p].back()
                S_sub[k - 1] = -A_ei_right[k] * all_w_l[p].back();
            }
        }
        schurSolver->matA.fillMat(S_diag, S_sub, S_sup);

    } else if (num_interfaces == 1) { // Caso P = 2 (n_interfaces = 1)
        // S è uno scalare 1x1
        schurScalarS = S_diag_ee[0]; // A_ee(0,0)
        // Contributo da dominio 0 (sinistra): A_ei_right[0] * w_r[0].back()
        schurScalarS -= A_ei_right[0] * all_w_r[0].back();
        // Contributo da dominio 1 (destra): A_ei_left[0] * w_l[1].front()
        schurScalarS -= A_ei_left[0] * all_w_l[1].front();
    }
}


// --- SOLVE (MODIFICATO) ---
std::vector<double> SchurSequentialSolver::solve(const std::vector<double>& f_global) {
    if (f_global.size() != nGlobal) {
        throw std::runtime_error("Dimensione f_global non corrisponde a nGlobal.");
    }

    std::vector<double> x_global(nGlobal);

    // Caso P = 1
    if (num_domains == 1) {
        localSolvers[0].setRhs(f_global);
        localSolvers[0].ThomaSolver();
        return localSolvers[0].getSolution();
    }

    // --- Casi P > 1 ---
    
    // Fasi 1, 2, 3 (invariate)
    // Fase 1: Partiziona f
    std::vector<std::vector<double>> f_i_all(num_domains);
    std::vector<double> f_e = getInterfaceVector(f_global);
    for (int p = 0; p < num_domains; ++p) {
        f_i_all[p] = getGlobalSubVector(f_global, domain_starts[p], domain_sizes[p]);
    }
    // Fase 2: Solve Locale 1 (z_p)
    std::vector<std::vector<double>> z_p_all(num_domains);
    for (int p = 0; p < num_domains; ++p) {
        localSolvers[p].setRhs(f_i_all[p]);
        localSolvers[p].ThomaSolver();
        z_p_all[p] = localSolvers[p].getSolution();
    }
    // Fase 3: Costruisci RHS di Schur (f*)
    std::vector<double> f_star = f_e; 
    for (int p = 0; p < num_domains; ++p) {
        if (p > 0) f_star[p-1] -= A_ei_left[p-1] * z_p_all[p][0];
        if (p < num_interfaces) f_star[p] -= A_ei_right[p] * z_p_all[p].back();
    }

    // --- Fase 4: Solve Interfaccia (MODIFICATA) ---
    std::vector<double> u_e;
    if (schurSolver) { // Caso P >= 3 (n_interfaces >= 2)
        schurSolver->setRhs(f_star);
        schurSolver->ThomaSolver();
        u_e = schurSolver->getSolution();
    } else if (num_interfaces == 1) { // Caso P = 2 (n_interfaces = 1)
        u_e.resize(1);
        if (schurScalarS == 0.0) {
            throw std::runtime_error("Divisione per zero: Schur scalar S è zero.");
        }
        u_e[0] = f_star[0] / schurScalarS;
    }
    // (Se P=1, non arriviamo mai qui)

    setInterfaceVector(x_global, u_e);

    // --- Fase 5: Solve Locale 2 (invariata) ---
    for (int p = 0; p < num_domains; ++p) {
        std::vector<double> g_p = f_i_all[p]; 
        if (p > 0) g_p[0] -= A_ie_left[p] * u_e[p-1];
        if (p < num_interfaces) g_p.back() -= A_ie_right[p] * u_e[p];

        localSolvers[p].setRhs(g_p);
        localSolvers[p].ThomaSolver();
        std::vector<double> u_i = localSolvers[p].getSolution();
        setGlobalSubVector(x_global, u_i, domain_starts[p]);
    }

    return x_global;
}


// --- Funzioni Helper (invariate) ---
// ... (tutte le funzioni helper restano identiche) ...
std::vector<double> SchurSequentialSolver::getGlobalSubVector(
    const std::vector<double>& globalVec, int start, int size) const 
{
    return std::vector<double>(globalVec.begin() + start, globalVec.begin() + start + size);
}
void SchurSequentialSolver::setGlobalSubVector(
    std::vector<double>& globalVec, const std::vector<double>& subVec, int start) const 
{
    std::copy(subVec.begin(), subVec.end(), globalVec.begin() + start);
}
std::vector<double> SchurSequentialSolver::getInterfaceVector(
    const std::vector<double>& globalVec) const 
{
    std::vector<double> interfaceVec(num_interfaces);
    for (int i = 0; i < num_interfaces; ++i) {
        interfaceVec[i] = globalVec[interface_indices[i]];
    }
    return interfaceVec;
}
void SchurSequentialSolver::setInterfaceVector(
    std::vector<double>& globalVec, const std::vector<double>& interfaceVec) const 
{
    for (int i = 0; i < num_interfaces; ++i) {
        globalVec[interface_indices[i]] = interfaceVec[i];
    }
}
std::vector<double> SchurSequentialSolver::createSparseVector(
    int size, int index, double value) const
{
    std::vector<double> vec(size, 0.0);
    if (index >= 0 && index < size) {
        vec[index] = value;
    }
    return vec;
}