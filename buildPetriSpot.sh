#!/bin/bash
# Build PetriSpot and its only dependency (a static libexpat) with CMake.
# Result: stripped petri32/64/128 (and kersconv) in website/.
#
# Optional argument on macOS: the Homebrew GCC to use, e.g. "gcc@13" (default).
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PREFIX="$SCRIPT_DIR/usr/local"
EXPAT_VERSION=2.6.4
JOBS=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)
EXTRA_CMAKE_FLAGS=""

# On macOS, brew install the requested GCC and use it (not clang).
if [ "$(uname)" = "Darwin" ]; then
  GCC_ARG="${1:-gcc@13}"
  GCC_VER="${GCC_ARG#gcc@}"
  brew install "gcc@${GCC_VER}" cmake
  export CC="gcc-${GCC_VER}"
  export CXX="g++-${GCC_VER}"
  export AR="gcc-ar-${GCC_VER}"
  export NM="gcc-nm-${GCC_VER}"
  export RANLIB="gcc-ranlib-${GCC_VER}"
fi

mkdir -p "$PREFIX"

# --- libexpat (static) ---
echo "=== Building libexpat $EXPAT_VERSION ==="
cd "$SCRIPT_DIR"
EXPAT_TAG="R_${EXPAT_VERSION//./_}"
if [ ! -f "expat-$EXPAT_VERSION.tar.gz" ]; then
  wget --progress=dot:mega "https://github.com/libexpat/libexpat/releases/download/$EXPAT_TAG/expat-$EXPAT_VERSION.tar.gz"
fi
rm -rf "expat-$EXPAT_VERSION"
tar xzf "expat-$EXPAT_VERSION.tar.gz"
cmake -S "expat-$EXPAT_VERSION" -B "expat-$EXPAT_VERSION/build" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DEXPAT_SHARED_LIBS=OFF -DEXPAT_BUILD_TOOLS=OFF -DEXPAT_BUILD_EXAMPLES=OFF \
    -DEXPAT_BUILD_TESTS=OFF -DEXPAT_BUILD_DOCS=OFF -DEXPAT_BUILD_PKGCONFIG=OFF
cmake --build "expat-$EXPAT_VERSION/build" -j"$JOBS"
cmake --install "expat-$EXPAT_VERSION/build"
rm -rf "expat-$EXPAT_VERSION" "expat-$EXPAT_VERSION.tar.gz"

# --- PetriSpot ---
echo "=== Building PetriSpot ==="
cd "$SCRIPT_DIR"
rm -rf build
cmake -S Petri -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$PREFIX" $EXTRA_CMAKE_FLAGS
cmake --build build -j"$JOBS"

# --- Package binaries ---
echo "=== Packaging binaries ==="
mkdir -p website
cp build/petri32 build/petri64 build/petri128 build/kersconv website/
strip website/petri32 website/petri64 website/petri128 website/kersconv

echo "=== Done. Binaries in $SCRIPT_DIR/website/ ==="
