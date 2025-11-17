#pragma once

#include <vector>
#include <memory>
#include "core/TridiagMat.hpp"
#include "numerics/LinearSys.hpp"

class SchurSequentialSolver {

public:
    SchurSequentialSolver(int globalSize, int numDomains, BoundaryType type);

    /**
     * @brief Fase 0: Pre-processing.
     * Costruisce le sottomatrici A_ii, estrae i coefficienti di interfaccia
     * e assembla la matrice di Schur 'S'.
     * @param A_global La matrice tridiagonale GLOBALE (N x N).
     */
    void PreProcess(const TridiagMat& A_global);

    /**
     * @brief Fase 1: Solve.
     * Risolve il sistema A*x = f usando l'algoritmo di Schur.
     * 'PreProcess' deve essere stato chiamato prima.
     * @param f_global Il vettore RHS globale (lunghezza N).
     * @return Il vettore soluzione globale x (lunghezza N).
     */
    std::vector<double> solve(const std::vector<double>& f_global);


private:
    // --- Variabili di Decomposizione ---
    int nGlobal;           // Dimensione totale N
    int num_domains;       // Numero di domini P
    int num_interfaces;    // Numero di interfacce (P-1)
    BoundaryType bType;

    std::vector<int> domain_sizes;  
    std::vector<int> domain_starts; 
    std::vector<int> interface_indices; 

    // --- Sistemi Locali e di Schur ---
    std::vector<LinearSys> localSolvers; // P solutori (per A_ii)

    // schurSolver viene creato SOLO se P >= 3 (num_interfaces >= 2)
    std::unique_ptr<LinearSys> schurSolver; // <-- MODIFICATO

    // schurScalarS viene usato SOLO se P = 2 (num_interfaces = 1)
    double schurScalarS = 0.0;

    // --- Dati delle Matrici di Interfaccia ---
    // A_ie: collega u_i (interno) a u_e (interfaccia)
    std::vector<double> A_ie_left;  // [P] coeff. 'a' (subdiag) che connettono a sinistra
    std::vector<double> A_ie_right; // [P] coeff. 'c' (supdiag) che connettono a destra

    // A_ei: collega u_e (interfaccia) a u_i (interno)
    std::vector<double> A_ei_left;  // [P-1] coeff. 'a' (subdiag) che connettono a sinistra
    std::vector<double> A_ei_right; // [P-1] coeff. 'c' (supdiag) che connettono a destra

    // --- Funzioni Helper Private ---
    std::vector<double> getGlobalSubVector(const std::vector<double>& globalVec, 
                                           int start, int size) const;
    void setGlobalSubVector(std::vector<double>& globalVec, 
                            const std::vector<double>& subVec, 
                            int start) const;
    std::vector<double> getInterfaceVector(const std::vector<double>& globalVec) const;
    void setInterfaceVector(std::vector<double>& globalVec, 
                            const std::vector<double>& interfaceVec) const;
    std::vector<double> createSparseVector(int size, int index, double value) const;
};