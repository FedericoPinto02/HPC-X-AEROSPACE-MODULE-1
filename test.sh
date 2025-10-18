set -euo pipefail

BUILD_DIR="build"
mkdir -p "$BUILD_DIR"
pushd "$BUILD_DIR" > /dev/null
CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=release"

cmake .. ${CMAKE_FLAGS}
make -j"$(nproc)"
popd > /dev/null

 "${BUILD_DIR}/testTridiag"
 "${BUILD_DIR}/test_derivatives"
 "${BUILD_DIR}/testFields"