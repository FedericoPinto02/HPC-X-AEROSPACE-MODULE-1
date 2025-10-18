#include "core/Fields.hpp"
#include <cmath>
#include <gtest/gtest.h>

// === FIELD CLASS ===

// --- getGrid test ---

TEST(FieldsTest, GetGrid_nonNullPointer) {
  std::shared_ptr<const Grid> gridPtr =
      std::make_shared<const Grid>(300, 100, 100, 0.3, 0.1, 0.1);
  Field field(gridPtr);

  auto returnedPtr = field.getGrid();

  EXPECT_NE(returnedPtr.get(), nullptr);
}

TEST(FieldsTest, GetGrid_correctPointer) {
  std::shared_ptr<const Grid> gridPtr =
      std::make_shared<const Grid>(300, 100, 100, 0.3, 0.1, 0.1);
  Field field(gridPtr);

  auto returnedPtr = field.getGrid();

  EXPECT_EQ(returnedPtr.get(), gridPtr.get());
}

TEST(FieldsTest, GetGrid_nullPointer) {
  Field field(nullptr);

  auto returnedPtr = field.getGrid();

  EXPECT_EQ(returnedPtr.get(), nullptr);
}