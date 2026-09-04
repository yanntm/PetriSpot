#ifndef SPARSE_MATRIX_IO_H
#define SPARSE_MATRIX_IO_H

#include "core/MatrixCol.h"
#include <bit>        // std::endian, std::byteswap (C++23)
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

/**
 * Binary format for sparse integer matrices (KERS format), see KERS.md.
 *
 * Header (16 bytes):
 *   magic   : 4 bytes  "KERS"
 *   version : 1 byte   (= 1)
 *   flags   : 1 byte   (reserved, = 0)
 *   nrows   : 4 bytes  uint32 LE
 *   ncols   : 4 bytes  uint32 LE
 *   padding : 2 bytes  (zero)
 *
 * Body: non-empty columns only, in ascending column order:
 *   col_idx    : 4 bytes  uint32 LE
 *   nnz        : 4 bytes  uint32 LE
 *   row_indices: nnz x 4 bytes uint32 LE (contiguous, sorted ascending)
 *   values     : nnz x 8 bytes int64  LE (contiguous)
 *
 * Terminator:
 *   0xFFFFFFFF : 4 bytes
 *
 * The stream overloads read or write exactly one such block and leave the
 * stream positioned after it, so blocks can be concatenated in a container
 * (see PNETIO.h). The file overloads keep the historical behaviour: errors
 * are reported on stderr and an empty matrix is returned.
 */

static const char     KERS_MAGIC[4] = {'K', 'E', 'R', 'S'};
static const uint32_t KERS_END      = 0xFFFFFFFF;

template<typename T>
class SparseMatrixIO {
public:

    /** Write one KERS block on os; false when the stream fails. */
    static bool write(const MatrixCol<T> &matrix, std::ostream &os) {
        uint32_t nrows = static_cast<uint32_t>(matrix.getRowCount());
        uint32_t ncols = static_cast<uint32_t>(matrix.getColumnCount());

        os.write(KERS_MAGIC, 4);
        os.put(1);              // version
        os.put(0);              // flags
        writeLE(os, nrows);
        writeLE(os, ncols);
        os.put(0); os.put(0);   // padding
        if (!os) return false;

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
            for (size_t i = 0; i < col.size(); ++i)
                appendLE(buf, static_cast<uint32_t>(col.keyAt(i)));
            for (size_t i = 0; i < col.size(); ++i)
                appendLE(buf, static_cast<int64_t>(col.valueAt(i)));
            os.write(buf.data(), static_cast<std::streamsize>(buf.size()));
            if (!os) return false;
        }
        writeLE(os, KERS_END);
        return static_cast<bool>(os);
    }

    /** Write a KERS file; errors on stderr, false on failure. */
    static bool write(const MatrixCol<T> &matrix, const std::string &filename) {
        std::ofstream ofs(filename, std::ios::binary | std::ios::out | std::ios::trunc);
        if (!ofs) {
            std::cerr << "Error: cannot open '" << filename << "' for writing\n";
            return false;
        }
        if (!write(matrix, ofs)) {
            std::cerr << "Error: write failed on '" << filename << "'\n";
            return false;
        }
        return true;
    }

    /**
     * Read one KERS block from is. Throws std::string on a malformed or
     * truncated block; label names the source in the message.
     */
    static MatrixCol<T> read(std::istream &is, uint32_t &nrows_out, uint32_t &ncols_out,
                             const std::string &label = "stream") {
        nrows_out = ncols_out = 0;
        uint8_t header[16];
        if (!is.read(reinterpret_cast<char *>(header), 16))
            throw std::string("Truncated KERS header in " + label);
        if (std::memcmp(header, KERS_MAGIC, 4) != 0)
            throw std::string("Bad KERS magic in " + label);
        uint8_t version = header[4];
        if (version != 1)
            throw std::string("Unsupported KERS version " + std::to_string(version) + " in " + label);

        uint32_t nrows, ncols;
        std::memcpy(&nrows, header + 6,  4);  nrows = toLittleEndian(nrows);
        std::memcpy(&ncols, header + 10, 4);  ncols = toLittleEndian(ncols);
        nrows_out = nrows;
        ncols_out = ncols;

        MatrixCol<T> matrix(nrows, ncols);
        std::vector<uint32_t> rows;
        std::vector<int64_t>  vals;
        uint32_t col_idx = KERS_END;
        for (;;) {
            if (!readLE(is, col_idx)) throw std::string("Missing KERS terminator in " + label);
            if (col_idx == KERS_END) break;
            if (col_idx >= ncols) throw std::string("Column index out of range in " + label);
            uint32_t nnz;
            if (!readLE(is, nnz)) throw std::string("Truncated nnz in " + label);
            SparseArray<T> &col = matrix.getColumn(col_idx);
            rows.resize(nnz);
            vals.resize(nnz);
            for (uint32_t i = 0; i < nnz; ++i)
                if (!readLE(is, rows[i])) throw std::string("Truncated row data in " + label);
            for (uint32_t i = 0; i < nnz; ++i)
                if (!readLE(is, vals[i])) throw std::string("Truncated value data in " + label);
            for (uint32_t i = 0; i < nnz; ++i) {
                if (rows[i] >= nrows) throw std::string("Row index out of range in " + label);
                checkOverflow(vals[i], col_idx, rows[i]);
                col.append(rows[i], static_cast<T>(vals[i]));
            }
        }
        return matrix;
    }

    /** Read a KERS file; errors on stderr and an empty matrix on failure. */
    static MatrixCol<T> read(const std::string &filename) {
        uint32_t nr, nc;
        return read(filename, nr, nc);
    }

    /** Read a KERS file and report its header dimensions. */
    static MatrixCol<T> read(const std::string &filename,
                             uint32_t &nrows_out, uint32_t &ncols_out) {
        nrows_out = ncols_out = 0;
        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs) {
            std::cerr << "Error: cannot open '" << filename << "' for reading\n";
            return MatrixCol<T>();
        }
        try {
            return read(ifs, nrows_out, ncols_out, "'" + filename + "'");
        } catch (const std::string &err) {
            std::cerr << "Error: " << err << "\n";
            nrows_out = ncols_out = 0;
            return MatrixCol<T>();
        }
    }

    // Zero-cost little-endian conversion (C++23)
    template<typename U>
    static U toLittleEndian(U v) {
        if constexpr (std::endian::native == std::endian::little)
            return v;
        else
            return std::byteswap(v);
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

private:

    template<typename U>
    static void appendLE(std::vector<char> &buf, U value) {
        U le = toLittleEndian(value);
        auto old_size = buf.size();
        buf.resize(old_size + sizeof(U));
        std::memcpy(buf.data() + old_size, &le, sizeof(U));
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
