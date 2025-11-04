#include "core/Fields.hpp"
#include <cmath>
#include <gtest/gtest.h>

// === FIELD CLASS ===

// Test Fixture
class FieldTestFixture : public ::testing::Test {
protected:
  const size_t Nx = 2, Ny = 3, Nz = 4;
  const size_t vectorSize = Nx * Ny * Nz;
  std::shared_ptr<const Grid> gridPtr;
  
  std::vector<Field::Scalar> m_v;
  Field testField;

  void SetUp() override {
    gridPtr = std::make_shared<const Grid>(Nx, Ny, Nz, 0.1, 0.1, 0.1);
    m_v.resize(vectorSize);
    for (size_t i = 0; i < vectorSize; i++)

      m_v.at(i) = (Field::Scalar)i;

    std::vector<Field::Scalar> setupDataCopy = m_v;
    testField.setup(gridPtr, setupDataCopy);
  }

  Field createOtherField(Field::Scalar initialValue = 0.0) {
    std::vector<Field::Scalar> setupData(vectorSize);
    
    std::fill(setupData.begin(), setupData.end(), initialValue);
    Field other;

    other.setup(gridPtr, setupData);
    return other;
  }

  Field createMismatchedField(size_t mismatchNx, size_t mismatchNy,
                              size_t mismatchNz) {
    size_t wrongSize = mismatchNx * mismatchNy * mismatchNz;
    std::shared_ptr<const Grid> wrongGridPtr = std::make_shared<const Grid>(
        mismatchNx, mismatchNy, mismatchNz, 0.1, 0.1, 0.1);
    std::vector<Field::Scalar> setupData(wrongSize);
    std::fill(setupData.begin(), setupData.end(), 1.0);

    Field other;
    other.setup(wrongGridPtr, setupData);
    return other;
  }
};

// === getGrid tests ===

TEST_F(FieldTestFixture, GetGrid_nonNullPointer) {
  EXPECT_NE(testField.getGrid().get(), nullptr);
}

TEST_F(FieldTestFixture, GetGrid_correctPointer) {
  EXPECT_EQ(testField.getGrid().get(), gridPtr.get());
}

// === access operator tests ===

TEST_F(FieldTestFixture, AccessOperator_firstElement) {
  EXPECT_DOUBLE_EQ(testField(0, 0, 0), m_v.at(0));
}

TEST_F(FieldTestFixture, AccessOperator_randomElementOne) {
  EXPECT_DOUBLE_EQ(testField(1, 1, 2), m_v.at(15));
}

TEST_F(FieldTestFixture, AccessOperator_randomElementTwo) {
  EXPECT_DOUBLE_EQ(testField(0, 2, 0), m_v.at(4));
}

TEST_F(FieldTestFixture, AccessOperator_lastElement) {
  EXPECT_DOUBLE_EQ(testField(testField.getGrid()->Nx - 1,
                             testField.getGrid()->Ny - 1,
                             testField.getGrid()->Nz - 1),
                   m_v.at(vectorSize - 1));
}

// === setup tests ===

TEST(FieldTest, Setup_gridPointerNull) {
  std::vector<Field::Scalar> m_v;
  Field field;

  EXPECT_THROW(field.setup(nullptr, m_v), std::invalid_argument);
}


TEST(FieldTest, Setup_vectorSizeMismatch) {
  std::shared_ptr<const Grid> gridPtr =
      std::make_shared<const Grid>(2, 3, 4, 0.1, 0.1, 0.1);
  size_t size_m_v = 2 * 3 * 4;
  std::vector<Field::Scalar> m_v;
  m_v.resize(size_m_v + 1);
  for (size_t i = 0; i < size_m_v; i++)
    m_v.at(i) = (Field::Scalar)i;
  Field field;

  EXPECT_THROW(field.setup(gridPtr, m_v), std::invalid_argument);
}


TEST_F(FieldTestFixture, Setup_correctGridSetup) {
  EXPECT_EQ(gridPtr.get(), testField.getGrid().get());
}

TEST_F(FieldTestFixture, Setup_correctVectorSetup) {
  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testField(i, j, k), m_v.at(index));
      }
}

// === reset test ===


TEST_F(FieldTestFixture, Reset_correctExecution) {
  testField.reset(1.0);
  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testField(i, j, k), 1.0);
      }
}

// === update tests ===
TEST_F(FieldTestFixture, Update_invalidArgumentNewVectorSize) {
  std::vector<Field::Scalar> m_v_new;
  m_v_new.resize(vectorSize - 2);
  std::fill(m_v_new.begin(), m_v_new.end(), 1.0);

  EXPECT_THROW(testField.update(m_v_new), std::invalid_argument);
}


TEST_F(FieldTestFixture, Update_correctExecution) {
  std::vector<Field::Scalar> m_v_new;
  m_v_new.resize(vectorSize);
  std::fill(m_v_new.begin(), m_v_new.end(), 2.5);
  testField.update(m_v_new);

  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
    
      for (size_t k = 0; k < Nz; k++)
        EXPECT_DOUBLE_EQ(testField(i, j, k), 2.5);
}

// === add (element wise) test ===

TEST_F(FieldTestFixture, AddElement_correctExecution) {
  Field::Scalar value = 2.0;
  testField.add(value);
  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testField(i, j, k), m_v.at(index) + value);
      }
}

// === add (field wise) tests ===

TEST_F(FieldTestFixture, Add_incorrectSizeGrid_x) {
  Field field_new = createMismatchedField(3, 3, 4);

  EXPECT_THROW(testField.add(field_new), std::invalid_argument);
}

TEST_F(FieldTestFixture, Add_incorrectSizeGrid_y) {
  Field field_new = createMismatchedField(2, 5, 4);

  EXPECT_THROW(testField.add(field_new), std::invalid_argument);
}

TEST_F(FieldTestFixture, Add_incorrectSizeGrid_z) {
  Field field_new = createMismatchedField(2, 3, 5);

  EXPECT_THROW(testField.add(field_new), std::invalid_argument);
}

TEST_F(FieldTestFixture, Add_correctExecution) {
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

TEST_F(FieldTestFixture, Multiply_correctExecution) {
  testField.multiply(0.5);
  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testField(i, j, k), m_v.at(index) * 0.5);
      }
}

// === VECTOR_FIELD CLASS ===

// Test Fixture
class VectorFieldTestFixture : public ::testing::Test {
protected:
  const size_t Nx = 2, Ny = 3, Nz = 4;
  const size_t vectorSize = Nx * Ny * Nz;
  std::shared_ptr<const Grid> gridPtr;
  std::vector<Field::Scalar> values_x, values_y, values_z;
  VectorField testVectorField;

  void SetUp() override {
    gridPtr = std::make_shared<const Grid>(Nx, Ny, Nz, 0.1, 0.1, 0.1);
    values_x.resize(vectorSize);
    values_y.resize(vectorSize);
    values_z.resize(vectorSize);
    for (size_t i = 0; i < vectorSize; i++) {
      values_x.at(i) = (Field::Scalar)i;
      values_y.at(i) = (Field::Scalar)i;
      values_z.at(i) = (Field::Scalar)i;
    }

    std::vector<Field::Scalar> setupDataCopy_x = values_x;
    std::vector<Field::Scalar> setupDataCopy_y = values_y;
    std::vector<Field::Scalar> setupDataCopy_z = values_z;
    testVectorField.setup(gridPtr, setupDataCopy_x, setupDataCopy_y,
                          setupDataCopy_z);
  }

  VectorField createOtherVectorField(Field::Scalar initialValue_x = 0.0,
                                     Field::Scalar initialValue_y = 1.0,
                                     Field::Scalar initialValue_z = 2.0) {
    std::vector<Field::Scalar> setupData_x(vectorSize);
    std::vector<Field::Scalar> setupData_y(vectorSize);
    std::vector<Field::Scalar> setupData_z(vectorSize);
    std::fill(setupData_x.begin(), setupData_x.end(), initialValue_x);
    std::fill(setupData_y.begin(), setupData_y.end(), initialValue_y);
    std::fill(setupData_z.begin(), setupData_z.end(), initialValue_z);
    VectorField other;

    other.setup(gridPtr, setupData_x, setupData_y, setupData_z);
    return other;
  }
};

// === getGrid tests ===

TEST_F(VectorFieldTestFixture, GetGrid_nonNullPointer) {
  EXPECT_NE(testVectorField.getGrid().get(), nullptr);
}

TEST_F(VectorFieldTestFixture, GetGrid_correctExecution) {
  EXPECT_EQ(testVectorField.getGrid(), gridPtr);
}

// === access operator vector x test ===

TEST_F(VectorFieldTestFixture, AccessOperatorX_correctExecution) {
  Field m_x = testVectorField.x();

  // to review
  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(m_x(i, j, k), values_x.at(index));
      }

  EXPECT_EQ(m_x.getGrid(), gridPtr);
}

// === access operator element on x vector tests ===

TEST_F(VectorFieldTestFixture, AccessOperatorElementOnX_firstElement) {
  size_t i = 0;
  size_t j = 0;
  size_t k = 0;
  Field::Scalar elem_x = testVectorField.x(i, j, k);

  size_t index = i + Nx * (j + Ny * k);
  EXPECT_DOUBLE_EQ(elem_x, values_x.at(index));
  EXPECT_EQ(testVectorField.x().getGrid(), gridPtr);
}

TEST_F(VectorFieldTestFixture, AccessOperatorElementOnX_randomElement) {
  size_t i = 1;
  size_t j = 0;
  size_t k = 3;
  Field::Scalar elem_x = testVectorField.x(i, j, k);

  size_t index = i + Nx * (j + Ny * k);
  EXPECT_DOUBLE_EQ(elem_x, values_x.at(index));
  EXPECT_EQ(testVectorField.x().getGrid(), gridPtr);
}

TEST_F(VectorFieldTestFixture, AccessOperatorElementOnX_lastElement) {
  size_t i = Nx - 1;
  size_t j = Ny - 1;
  size_t k = Nz - 1;
  Field::Scalar elem_x = testVectorField.x(i, j, k);

  size_t index = i + Nx * (j + Ny * k);
  EXPECT_DOUBLE_EQ(elem_x, values_x.at(index));
  EXPECT_EQ(testVectorField.x().getGrid(), gridPtr);
}

// === access operator vector y test ===

TEST_F(VectorFieldTestFixture, AccessOperatorY_correctExecution) {
  Field m_y = testVectorField.y();

  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(m_y(i, j, k), values_y.at(index));
      }

  EXPECT_EQ(m_y.getGrid(), gridPtr);
}

// === access operator element on y vector tests ===

TEST_F(VectorFieldTestFixture, AccessOperatorElementOnY_firstElement) {
  size_t i = 0;
  size_t j = 0;
  size_t k = 0;
  Field::Scalar elem_y = testVectorField.y(i, j, k);

  size_t index = i + Nx * (j + Ny * k);
  EXPECT_DOUBLE_EQ(elem_y, values_y.at(index));
  EXPECT_EQ(testVectorField.y().getGrid(), gridPtr);
}

TEST_F(VectorFieldTestFixture, AccessOperatorElementOnY_randomElement) {
  size_t i = 0;
  size_t j = 2;
  size_t k = 1;
  Field::Scalar elem_y = testVectorField.y(i, j, k);

  size_t index = i + Nx * (j + Ny * k);
  EXPECT_DOUBLE_EQ(elem_y, values_y.at(index));
  EXPECT_EQ(testVectorField.y().getGrid(), gridPtr);
}

TEST_F(VectorFieldTestFixture, AccessOperatorElementOnY_lastElement) {
  size_t i = Nx - 1;
  size_t j = Ny - 1;
  size_t k = Nz - 1;
  Field::Scalar elem_y = testVectorField.y(i, j, k);

  size_t index = i + Nx * (j + Ny * k);
  EXPECT_DOUBLE_EQ(elem_y, values_y.at(index));
  EXPECT_EQ(testVectorField.y().getGrid(), gridPtr);
}

// === access operator vector z test ===

TEST_F(VectorFieldTestFixture, AccessOperatorZ_correctExecution) {
  Field m_z = testVectorField.z();

  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(m_z(i, j, k), values_z.at(index));
      }

  EXPECT_EQ(m_z.getGrid(), gridPtr);
}

// === access operator element on z vector tests ===

TEST_F(VectorFieldTestFixture, AccessOperatorElementOnZ_firstElement) {
  size_t i = 0;
  size_t j = 0;
  size_t k = 0;
  Field::Scalar elem_z = testVectorField.z(i, j, k);

  size_t index = i + Nx * (j + Ny * k);
  EXPECT_DOUBLE_EQ(elem_z, values_z.at(index));
  EXPECT_EQ(testVectorField.z().getGrid(), gridPtr);
}

TEST_F(VectorFieldTestFixture, AccessOperatorElementOnZ_randomElement) {
  size_t i = 0;
  size_t j = 3;
  size_t k = 1;
  Field::Scalar elem_z = testVectorField.z(i, j, k);

  size_t index = i + Nx * (j + Ny * k);
  EXPECT_DOUBLE_EQ(elem_z, values_z.at(index));
  EXPECT_EQ(testVectorField.z().getGrid(), gridPtr);
}

TEST_F(VectorFieldTestFixture, AccessOperatorElementOnZ_lastElement) {
  size_t i = Nx - 1;
  size_t j = Ny - 1;
  size_t k = Nz - 1;
  Field::Scalar elem_z = testVectorField.z(i, j, k);

  size_t index = i + Nx * (j + Ny * k);
  EXPECT_DOUBLE_EQ(elem_z, values_z.at(index));
  EXPECT_EQ(testVectorField.z().getGrid(), gridPtr);
}

// === setup tests ===

TEST(VectorFieldTest, Setup_gridPointerNull) {
  std::vector<Field::Scalar> values_x, values_y, values_z;
  VectorField vectorField;

  EXPECT_THROW(vectorField.setup(nullptr, values_x, values_y, values_z),
               std::invalid_argument);
}

TEST(VectorFieldTest, Setup_correctExecution) {
  const size_t Nx = 2, Ny = 3, Nz = 4;
  const size_t vectorSize = Nx * Ny * Nz;
  std::shared_ptr<const Grid> gridPtr;
  std::vector<Field::Scalar> values_x, values_y, values_z;
  VectorField testVectorField;

  gridPtr = std::make_shared<const Grid>(Nx, Ny, Nz, 0.1, 0.1, 0.1);
  values_x.resize(vectorSize);
  values_y.resize(vectorSize);
  values_z.resize(vectorSize);
  for (size_t i = 0; i < vectorSize; i++) {
    values_x.at(i) = (Field::Scalar)i;
    values_y.at(i) = (Field::Scalar)i;
    values_z.at(i) = (Field::Scalar)i;
  }

  testVectorField.setup(gridPtr, values_x, values_y, values_z);

  EXPECT_EQ(testVectorField.getGrid(), gridPtr);
  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testVectorField.x()(i, j, k), values_x.at(index));
        EXPECT_DOUBLE_EQ(testVectorField.y()(i, j, k), values_y.at(index));
        EXPECT_DOUBLE_EQ(testVectorField.z()(i, j, k), values_z.at(index));
      }
  EXPECT_EQ(testVectorField.x().getGrid(), gridPtr);
  EXPECT_EQ(testVectorField.y().getGrid(), gridPtr);
  EXPECT_EQ(testVectorField.z().getGrid(), gridPtr);
}

// === update test ===

TEST_F(VectorFieldTestFixture, Update_correctExecution) {
  std::vector<Field::Scalar> values_x_new, values_y_new, values_z_new;
  values_x_new.resize(vectorSize);
  values_y_new.resize(vectorSize);
  values_z_new.resize(vectorSize);
  std::fill(values_x_new.begin(), values_x_new.end(), 0.5);
  std::fill(values_y_new.begin(), values_y_new.end(), 1.0);
  std::fill(values_z_new.begin(), values_z_new.end(), 1.5);

  testVectorField.update(values_x_new, values_y_new, values_z_new);

  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        EXPECT_DOUBLE_EQ(testVectorField.x()(i, j, k), 0.5);
        EXPECT_DOUBLE_EQ(testVectorField.y()(i, j, k), 1.0);
        EXPECT_DOUBLE_EQ(testVectorField.z()(i, j, k), 1.5);
      }
}

// === add (element wise) test ===

TEST_F(VectorFieldTestFixture, addElement_correctExecution) {
  Field::Scalar value = -2.5;
  testVectorField.add(value);

  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testVectorField.x()(i, j, k),
                         values_x.at(index) + value);
        EXPECT_DOUBLE_EQ(testVectorField.y()(i, j, k),
                         values_y.at(index) + value);
        EXPECT_DOUBLE_EQ(testVectorField.z()(i, j, k),
                         values_z.at(index) + value);
      }
}

// === add (field wise) test ===

TEST_F(VectorFieldTestFixture, addField_correctExecution) {
  VectorField otherVectorField = createOtherVectorField();
  testVectorField.add(otherVectorField);

  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testVectorField.x()(i, j, k),
                         values_x.at(index) + otherVectorField.x()(i, j, k));
        EXPECT_DOUBLE_EQ(testVectorField.y()(i, j, k),
                         values_y.at(index) + otherVectorField.y()(i, j, k));
        EXPECT_DOUBLE_EQ(testVectorField.z()(i, j, k),
                         values_z.at(index) + otherVectorField.z()(i, j, k));
      }
}

// === multiply test ===

TEST_F(VectorFieldTestFixture, multiply_correctExecution) {
  Field::Scalar value = 0.75;
  testVectorField.multiply(value);

  for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
      for (size_t k = 0; k < Nz; k++) {
        size_t index = i + Nx * (j + Ny * k);
        EXPECT_DOUBLE_EQ(testVectorField.x()(i, j, k),
                         values_x.at(index) * value);
        EXPECT_DOUBLE_EQ(testVectorField.y()(i, j, k),
                         values_y.at(index) * value);
        EXPECT_DOUBLE_EQ(testVectorField.z()(i, j, k),
                         values_z.at(index) * value);
      }
}