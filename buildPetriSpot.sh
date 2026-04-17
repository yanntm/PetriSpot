#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PREFIX="$SCRIPT_DIR/usr/local"

# On macOS, brew installs glibtoolize instead of libtoolize
if [ "$(uname)" = "Darwin" ]; then
  export LIBTOOLIZE=glibtoolize
fi

mkdir -p "$PREFIX"

# --- libExpat ---
echo "=== Building libExpat ==="
cd "$SCRIPT_DIR"
wget --progress=dot:mega https://github.com/libexpat/libexpat/archive/R_2_2_4.tar.gz
tar xzf R_2_2_4.tar.gz
cd libexpat-R_2_2_4/expat/
./buildconf.sh
./configure --prefix="$PREFIX" --without-xmlwf
make -j$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)
make install
cd "$SCRIPT_DIR"
rm -rf libexpat-R_2_2_4 R_2_2_4.tar.gz

# --- PetriSpot ---
echo "=== Building PetriSpot ==="
cd "$SCRIPT_DIR/Petri"
autoreconf -vfi
./configure \
    --prefix="$PREFIX" \
    --with-libexpat="$PREFIX" \
    CPPFLAGS="-I$PREFIX/include -DNDEBUG" \
    LDFLAGS="-L$PREFIX/lib" \
    || { cat config.log ; exit 1 ; }
make -j$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)

# --- Package binaries ---
echo "=== Packaging binaries ==="
cd "$SCRIPT_DIR"
mkdir -p website
cp Petri/src/petri32 Petri/src/petri64 Petri/src/petri128 website/
strip website/petri32 website/petri64 website/petri128

echo "=== Done. Binaries in $SCRIPT_DIR/website/ ==="
