#include <gtest/gtest.h>
#include <cmath>
#include "core/TridiagMat.hpp"

TEST(TridiagMatTest, GetElementInBounds) {
    // 1. Setup: Crea e popola una matrice 4x4
    TridiagMat matrix(4);
    matrix.fillMat(
        {10.0, 20.0, 30.0, 40.0},   // Diagonale
        {51.0, 61.0, 71.0},         // Sotto-diagonale
        {15.0, 25.0, 35.0}          // Sopra-diagonale
    );

    // 2. Asserzioni: Verifica i valori attesi
    
    // Controlla la diagonale principale
    EXPECT_DOUBLE_EQ(matrix.getElement(1, 1), 20.0);
    
}