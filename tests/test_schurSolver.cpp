#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <stdexcept>
#include <numeric> 

#include "numerics/LinearSys.hpp"
#include "numerics/SchurSequentialSolver.hpp"
#include "core/TridiagMat.hpp"

class SchurSolverTestFixture : public ::testing::TestWithParam<std::tuple<int, int>> {
protected:
    int N_global;
    int P_domains;

    std::unique_ptr<TridiagMat> A_global_ptr; 
    std::vector<double> f_global;
    BoundaryType bType = BoundaryType::Normal; 

    std::unique_ptr<LinearSys> directSolver;
    std::unique_ptr<SchurSequentialSolver> schurSolver;

    void createTestMatrix(TridiagMat& mat) {
        int n = mat.getSize();
        if (n == 0) return;

        std::vector<double> diag(n, 2.0);
        std::vector<double> sub(n - 1, -1.0);
        std::vector<double> sup(n - 1, -1.0);

        if (n > 0) {
            diag[0] = 1.0;
            if (n > 1) sup[0] = 0.0;
        }
        if (n > 1) {
            diag[n-1] = 1.0;
            sub[n-2] = 0.0;
        }

        mat.fillMat(diag, sub, sup);
    }

    void createTestRhs(std::vector<double>& rhs) {
        std::fill(rhs.begin(), rhs.end(), 0.0);
        if (rhs.size() > 2) {
            rhs[rhs.size() / 2] = 1.0;
        } else if (!rhs.empty()) {
            rhs[0] = 1.0;
        }
    }

    void SetUp() override {
        N_global = std::get<0>(GetParam());
        P_domains = std::get<1>(GetParam());
        
        // 1. Inizializza matrici e vettori globali
        A_global_ptr = std::make_unique<TridiagMat>(N_global);
        f_global.resize(N_global);

        createTestMatrix(*A_global_ptr);
        createTestRhs(f_global);

        // 2. Inizializza il solutore diretto
        directSolver = std::make_unique<LinearSys>(N_global, bType);
        
        // --- QUESTA È LA CORREZIONE ---
        // Errore: non possiamo fare l'assegnazione
        // directSolver->matA = *A_global_ptr; // DELETED!
        
        // Correzione: Usiamo fillMat() per copiare i contenuti
        // (Abbiamo ancora bisogno dell'accesso 'friend' a LinearSys 
        // per accedere al membro privato 'matA' qui)
        directSolver->matA.fillMat(
            A_global_ptr->getDiag(0),  // Copia la diagonale
            A_global_ptr->getDiag(-1), // Copia la sottodiagonale
            A_global_ptr->getDiag(1)   // Copia la sopradiagonale
        );
        // ---------------------------------

        // 3. Inizializza il solutore di Schur e fai il PreProcess
        schurSolver = std::make_unique<SchurSequentialSolver>(N_global, P_domains, bType);
        schurSolver->PreProcess(*A_global_ptr);
    }
};

// ... (Il resto del file rimane invariato) ...

TEST_P(SchurSolverTestFixture, SolutionAgreesWithDirectSolve) {
    // 1. Risolvi con il metodo diretto (seriale)
    directSolver->setRhs(f_global);
    directSolver->ThomaSolver();
    std::vector<double> x_direct = directSolver->getSolution();

    // 2. Risolvi con il metodo di Schur
    std::vector<double> x_schur = schurSolver->solve(f_global);

    // 3. Confronta
    ASSERT_EQ(x_direct.size(), x_schur.size());
    for (size_t i = 0; i < N_global; ++i) {
        EXPECT_NEAR(x_direct[i], x_schur[i], 1e-12)
            << "I vettori differiscono all'indice " << i;
    }
}

INSTANTIATE_TEST_SUITE_P(
    SchurVerification,
    SchurSolverTestFixture,
    ::testing::Values(
        std::make_tuple(50, 4), 
        std::make_tuple(100, 10),
        std::make_tuple(17, 5),  
        std::make_tuple(20, 1),  
        std::make_tuple(2, 1),   
        std::make_tuple(5, 2)    
    )
);

TEST(SchurSolverConstructorTest, HandlesInvalidInputs) {
    EXPECT_THROW(SchurSequentialSolver(10, 0, BoundaryType::Normal), std::invalid_argument);
    EXPECT_THROW(SchurSequentialSolver(5, 10, BoundaryType::Normal), std::invalid_argument);
    EXPECT_NO_THROW(SchurSequentialSolver(10, 1, BoundaryType::Normal));
}