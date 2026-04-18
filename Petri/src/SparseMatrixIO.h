#ifndef SPARSE_MATRIX_IO_H
#define SPARSE_MATRIX_IO_H

#include "MatrixCol.h"
#include <bit>        // std::endian, std::byteswap (C++23)
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

/**
 * Binary format for sparse integer matrices (KERS format).
 *
 * Header (16 bytes):
 *   magic   : 4 bytes  "KERS"
 *   version : 1 byte   (= 1)
 *   flags   : 1 byte   (reserved, = 0)
 *   nrows   : 4 bytes  uint32 LE
 *   ncols   : 4 bytes  uint32 LE
 *   padding : 2 bytes  (zero)
 *
 * Body — non-empty columns only, in ascending column order:
 *   col_idx    : 4 bytes  uint32 LE
 *   nnz        : 4 bytes  uint32 LE
 *   row_indices: nnz × 4 bytes uint32 LE (contiguous, sorted ascending)
 *   values     : nnz × 8 bytes int64  LE (contiguous)
 *
 * Terminator:
 *   0xFFFFFFFF : 4 bytes
 */

static const char     KERS_MAGIC[4] = {'K', 'E', 'R', 'S'};
static const uint32_t KERS_END      = 0xFFFFFFFF;

template<typename T>
class SparseMatrixIO {
public:

    static bool write(const MatrixCol<T> &matrix, const std::string &filename) {
        std::ofstream ofs(filename, std::ios::binary | std::ios::out | std::ios::trunc);
        if (!ofs) {
            std::cerr << "Error: cannot open '" << filename << "' for writing\n";
            return false;
        }

        uint32_t nrows = static_cast<uint32_t>(matrix.getRowCount());
        uint32_t ncols = static_cast<uint32_t>(matrix.getColumnCount());

        // Header (version 1)
        ofs.write(KERS_MAGIC, 4);
        ofs.put(1);              // version
        ofs.put(0);              // flags
        writeLE(ofs, nrows);
        writeLE(ofs, ncols);
        ofs.put(0); ofs.put(0);  // padding

        if (!ofs) {
            std::cerr << "Error: failed to write header to '" << filename << "'\n";
            return false;
        }

        // Pre-scan: size buffer to largest column to avoid repeated reallocations
        const auto &cols = matrix.getColumns();
        size_t maxColBytes = 0;
        for (uint32_t ci = 0; ci < ncols; ++ci) {
            size_t nnz = cols[ci].size();
            if (nnz == 0) continue;
            size_t bytes = 8 + nnz * (4 + 8);   // col_idx+nnz header + rows + values
            if (bytes > maxColBytes) maxColBytes = bytes;
        }

        std::vector<char> buf;
        if (maxColBytes > 0) buf.reserve(maxColBytes);

        // Body: one write per non-empty column
        for (uint32_t ci = 0; ci < ncols; ++ci) {
            const auto &col = cols[ci];
            if (col.size() == 0) continue;

            buf.clear();
            appendLE(buf, ci);
            appendLE(buf, static_cast<uint32_t>(col.size()));

            // Contiguous row indices
            for (size_t i = 0; i < col.size(); ++i)
                appendLE(buf, static_cast<uint32_t>(col.keyAt(i)));

            // Contiguous values
            for (size_t i = 0; i < col.size(); ++i)
                appendLE(buf, static_cast<int64_t>(col.valueAt(i)));

            ofs.write(buf.data(), static_cast<std::streamsize>(buf.size()));
            if (!ofs) {
                std::cerr << "Error: write failed for column " << ci << "\n";
                return false;
            }
        }

        writeLE(ofs, KERS_END);
        if (!ofs) {
            std::cerr << "Error: failed to write terminator\n";
            return false;
        }

        return true;  // RAII closes the file
    }

    // Read a KERS file (v1 or v2). Returns an empty matrix on error.
    static MatrixCol<T> read(const std::string &filename) {
        uint32_t nr, nc;
        return readImpl(filename, nr, nc);
    }

    // Read and also return header dimensions separately.
    static MatrixCol<T> read(const std::string &filename,
                             uint32_t &nrows_out, uint32_t &ncols_out) {
        return readImpl(filename, nrows_out, ncols_out);
    }

private:

    // Zero-cost little-endian conversion (C++23)
    template<typename U>
    static U toLittleEndian(U v) {
        if constexpr (std::endian::native == std::endian::little)
            return v;
        else
            return std::byteswap(v);
    }

    template<typename U>
    static void appendLE(std::vector<char> &buf, U value) {
        U le = toLittleEndian(value);
        auto old_size = buf.size();
        buf.resize(old_size + sizeof(U));
        std::memcpy(buf.data() + old_size, &le, sizeof(U));
    }

    static void writeLE(std::ostream &os, uint32_t value) {
        uint32_t le = toLittleEndian(value);
        os.write(reinterpret_cast<const char *>(&le), 4);
    }

    template<typename U>
    static bool readLE(std::istream &is, U &out) {
        U raw;
        if (!is.read(reinterpret_cast<char *>(&raw), sizeof(U))) return false;
        out = toLittleEndian(raw);
        return true;
    }

    static MatrixCol<T> readImpl(const std::string &filename,
                                 uint32_t &nrows_out, uint32_t &ncols_out) {
        nrows_out = ncols_out = 0;

        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs) {
            std::cerr << "Error: cannot open '" << filename << "' for reading\n";
            return MatrixCol<T>();
        }

        uint8_t header[16];
        if (!ifs.read(reinterpret_cast<char *>(header), 16)) {
            std::cerr << "Error: truncated header in '" << filename << "'\n";
            return MatrixCol<T>();
        }
        if (std::memcmp(header, KERS_MAGIC, 4) != 0) {
            std::cerr << "Error: bad magic in '" << filename << "'\n";
            return MatrixCol<T>();
        }

        uint8_t version = header[4];
        if (version != 1) {
            std::cerr << "Error: unsupported KERS version " << (int)version << "\n";
            return MatrixCol<T>();
        }

        uint32_t nrows, ncols;
        std::memcpy(&nrows, header + 6,  4);  nrows = toLittleEndian(nrows);
        std::memcpy(&ncols, header + 10, 4);  ncols = toLittleEndian(ncols);
        nrows_out = nrows;
        ncols_out = ncols;

        MatrixCol<T> matrix(nrows, ncols);

        std::vector<uint32_t> rows;
        std::vector<int64_t>  vals;

        uint32_t col_idx = KERS_END;
        while (readLE(ifs, col_idx)) {
            if (col_idx == KERS_END) break;

            uint32_t nnz;
            if (!readLE(ifs, nnz))
                throw "Truncated nnz in KERS file.";

            SparseArray<T> &col = matrix.getColumn(col_idx);

            // rows block, then values block
            rows.resize(nnz);
            vals.resize(nnz);
            for (uint32_t i = 0; i < nnz; ++i)
                if (!readLE(ifs, rows[i])) throw "Truncated row data in KERS file.";
            for (uint32_t i = 0; i < nnz; ++i)
                if (!readLE(ifs, vals[i])) throw "Truncated value data in KERS file.";
            for (uint32_t i = 0; i < nnz; ++i) {
                checkOverflow(vals[i], col_idx, rows[i]);
                col.append(rows[i], static_cast<T>(vals[i]));
            }
        }

        return matrix;
    }

    static void checkOverflow(int64_t val, uint32_t col, uint32_t row) {
        if (val > static_cast<int64_t>(std::numeric_limits<T>::max()) ||
            val < static_cast<int64_t>(std::numeric_limits<T>::min())) {
            std::cerr << "Warning: value " << val
                      << " at col=" << col << " row=" << row
                      << " overflows T (" << sizeof(T) << " bytes); truncated.\n";
        }
    }
};

#endif // SPARSE_MATRIX_IO_H
