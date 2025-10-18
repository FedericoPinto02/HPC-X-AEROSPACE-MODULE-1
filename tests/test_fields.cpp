#include "core/Fields.hpp"
#include <cmath>
#include <gtest/gtest.h>

// === FIELD CLASS ===

// Test Fixture
class FieldsTestFixture : public ::testing::Test {
protected:
  const size_t Nx = 2, Ny = 3, Nz = 4;
  const size_t vectorSize = Nx * Ny * Nz;
  std::shared_ptr<const Grid> gridPtr;
  std::vector<double> m_v;
  Field testField;

  void SetUp() override {
    gridPtr = std::make_shared<const Grid>(Nx, Ny, Nz, 0.1, 0.1, 0.1);
    m_v.resize(vectorSize);
    for (size_t i = 0; i < vectorSize; i++)
      m_v.at(i) = (double)i;

    std::vector<double> setupDataCopy = m_v;
    testField.setup(gridPtr, setupDataCopy);
  }

  Field createOtherField(double initialValue = 0.0) {
    std::vector<double> setupData(vectorSize);
    std::fill(setupData.begin(), setupData.end(), initialValue);
    Field other;

    other.setup(gridPtr, setupData);
    return other;
  }

  Field createMismatchedField(int mismatchNx, int mismatchNy, int mismatchNz) {
    size_t wrongSize = mismatchNx * mismatchNy * mismatchNz;
    std::shared_ptr<const Grid> wrongGridPtr = std::make_shared<const Grid>(
        mismatchNx, mismatchNy, mismatchNz, 0.1, 0.1, 0.1);
    std::vector<double> setupData(wrongSize);
    std::fill(setupData.begin(), setupData.end(), 1.0);

    Field other;
    other.setup(wrongGridPtr, setupData);
    return other;
  }
};

// === getGrid tests ===

TEST_F(FieldsTestFixture, GetGrid_nonNullPointer) {
  EXPECT_NE(testField.getGrid().get(), nullptr);
}

TEST_F(FieldsTestFixture, GetGrid_correctPointer) {
  EXPECT_EQ(testField.getGrid().get(), gridPtr.get());
}

// === access operator tests ===

TEST_F(FieldsTestFixture, AccessOperator_firstElement) {
  EXPECT_DOUBLE_EQ(testField(0, 0, 0), m_v.at(0));
}

TEST_F(FieldsTestFixture, AccessOperator_randomElementOne) {
  EXPECT_DOUBLE_EQ(testField(1, 1, 2), m_v.at(15));
}

TEST_F(FieldsTestFixture, AccessOperator_randomElementTwo) {
  EXPECT_DOUBLE_EQ(testField(0, 2, 0), m_v.at(4));
}

TEST_F(FieldsTestFixture, AccessOperator_lastElement) {
  EXPECT_DOUBLE_EQ(testField(testField.getGrid()->Nx - 1,
                             testField.getGrid()->Ny - 1,
                             testField.getGrid()->Nz - 1),
                   m_v.at(vectorSize - 1));
}

// === setup tests ===

TEST(FieldsTest, Setup_gridPointerNull) {
  std::vector<double> m_v;
  Field field;

  EXPECT_THROW(field.setup(nullptr, m_v), std::invalid_argument);
}

TEST(FieldsTest, Setup_vectorSizeMismatch) {
  std::shared_ptr<const Grid> gridPtr =
      std::make_shared<const Grid>(2, 3, 4, 0.1, 0.1, 0.1);
  size_t size_m_v = 2 * 3 * 4;
  std::vector<double> m_v;
  m_v.resize(size_m_v + 1);
  for (size_t i = 0; i < size_m_v; i++)
    m_v.at(i) = (double)i;
  Field field;

  EXPECT_THROW(field.setup(gridPtr, m_v), std::invalid_argument);
}

TEST_F(FieldsTestFixture, Setup_correctGridSetup) {
  EXPECT_EQ(gridPtr.get(), testField.getGrid().get());
}

TEST_F(FieldsTestFixture, Setup_correctVectorSetup) {
  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testField(i, j, k), m_v.at(index));
      }
}

// === reset test ===

TEST_F(FieldsTestFixture, Reset_correctExecution) {
  testField.reset(1.0);
  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testField(i, j, k), 1.0);
      }
}

// === update tests ===

TEST_F(FieldsTestFixture, Update_invalidArgumentNewVectorSize) {
  std::vector<double> m_v_new;
  m_v_new.resize(vectorSize - 2);
  std::fill(m_v_new.begin(), m_v_new.end(), 1.0);

  EXPECT_THROW(testField.update(m_v_new), std::invalid_argument);
}

TEST_F(FieldsTestFixture, Update_correctExecution) {
  std::vector<double> m_v_new;
  m_v_new.resize(vectorSize);
  std::fill(m_v_new.begin(), m_v_new.end(), 2.5);
  testField.update(m_v_new);

  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testField(i, j, k), 2.5);
      }
}

// === add (element wise) test ===

TEST_F(FieldsTestFixture, AddElement_correctExecution) {
  testField.add(2.0);
  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testField(i, j, k), m_v.at(index) + 2.0);
      }
}

// === add (field wise) tests ===

TEST_F(FieldsTestFixture, Add_incorrectSizeGrid_x) {
  Field field_new = createMismatchedField(3, 3, 4);

  EXPECT_THROW(testField.add(field_new), std::invalid_argument);
}

TEST_F(FieldsTestFixture, Add_incorrectSizeGrid_y) {
  Field field_new = createMismatchedField(2, 5, 4);

  EXPECT_THROW(testField.add(field_new), std::invalid_argument);
}

TEST_F(FieldsTestFixture, Add_incorrectSizeGrid_z) {
  Field field_new = createMismatchedField(2, 3, 5);

  EXPECT_THROW(testField.add(field_new), std::invalid_argument);
}

TEST_F(FieldsTestFixture, Add_correctExecution) {
  Field field_new = createOtherField(1.5);
  testField.add(field_new);

  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testField(i, j, k), m_v.at(index) + 1.5);
      }
}

// === multiply test ===

TEST_F(FieldsTestFixture, Multiply_correctExecution) {
  testField.multiply(0.5);
  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testField(i, j, k), m_v.at(index) * 0.5);
      }
}