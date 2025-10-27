#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath>

#include "core/Fields.hpp"
#include "numerics/LinearSys.hpp"
#include "core/Mesh.hpp"


// Fixture semplificata per test di accesso e operazioni scalari
class FieldSimpleTestFixture : public ::testing::Test {
protected:
    // --- Configurazione Griglia Principale ---
    const size_t Nx = 2, Ny = 3, Nz = 4;
    const size_t vectorSize = Nx * Ny * Nz; // 24
    std::shared_ptr<const Grid> gridPtr;

    // --- Campo di Test Principale ---
    Field testField; // Campo popolato con 0, 1, 2, ...

    // --- Campo per Test di Errore ---
    std::shared_ptr<const Grid> mismatchedGridPtr;
    Field mismatchedField; // Ha una griglia diversa

    /**
     * @brief Inizializza gli oggetti prima di ogni test.
     */
    void SetUp() override {
        // 1. Inizializza la griglia principale
        gridPtr = std::make_shared<const Grid>(Nx, Ny, Nz, 0.1, 0.1, 0.1);

        // 2. Inizializza testField (0, 1, 2, ..., 23)
        std::vector<Field::Scalar> initialData(vectorSize);
        for (size_t i = 0; i < vectorSize; i++) {
            initialData[i] = static_cast<Field::Scalar>(i);
        }
        testField.setup(gridPtr, initialData);

        // 3. Inizializza la griglia e il campo per i test di errore
        mismatchedGridPtr = std::make_shared<const Grid>(5, 5, 5, 0.1, 0.1, 0.1);
        std::vector<Field::Scalar> mismatchedData(125, 1.0);
        mismatchedField.setup(mismatchedGridPtr, mismatchedData);
    }
};

/**
 * @test SetupCorrectlyInitializesGrid
 * @brief Verifica che il campo sia stato inizializzato correttamente
 * con la griglia e la dimensione corrette definite in SetUp().
 *
 * Questo è un test "ovvio" perché non controlla nessun valore di
 * elemento, ma solo che il costruttore/setup del campo abbia
 * memorizzato correttamente la sua griglia e calcolato la sua dimensione totale.
 */
TEST_F(FieldSimpleTestFixture, SetupCorrectlyInitializesGrid)
{
    // 1. Verifica che il campo abbia effettivamente una griglia
    ASSERT_NE(testField.getGrid(), nullptr);

    // 2. Verifica che la griglia memorizzata abbia le dimensioni corrette
    //    (Nx, Ny, Nz sono ereditati dalla fixture)
    EXPECT_EQ(testField.getGrid()->Nx, Nx);
    EXPECT_EQ(testField.getGrid()->Ny, Ny);
    EXPECT_EQ(testField.getGrid()->Nz, Nz);
    
    // 3. Verifica che il campo 'sappia' la sua dimensione totale
    //    (vectorSize è ereditato dalla fixture)
    EXPECT_EQ(testField.getGrid()->Nx, vectorSize);
}