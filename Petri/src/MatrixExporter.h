#ifndef MATRIX_EXPORTER_H
#define MATRIX_EXPORTER_H

#include "MatrixCol.h"
#include <fstream>
#include <string>

template<typename T>
class MatrixExporter {
public:
    static void exportMatrix(const MatrixCol<T>& matrix, const std::string& filename) {
        std::ofstream out(filename);
        if (!out) {
            std::cerr << "Error: Could not open " << filename << " for writing" << std::endl;
            return;
        }

        // Write header: nbRow nbCol
        size_t rows = matrix.getRowCount();
        size_t cols = matrix.getColumnCount();
        out << rows << " " << cols << "\n";

        // Transpose for row-wise output
        MatrixCol<T> transposed = matrix.transpose();

        // Iterate over columns of transposed matrix (original rows)
        const auto& columns = transposed.getColumns();
        for (size_t col = 0; col < transposed.getColumnCount(); ++col) {  // col is original row index
            const SparseArray<T>& sparse_col = columns[col];
            if (sparse_col.size() > 0) {  // Only write non-empty rows
                out << col;  // Original row index
                for (size_t i = 0; i < sparse_col.size(); ++i) {
                    size_t orig_col = sparse_col.keyAt(i);  // Original column index
                    T val = sparse_col.valueAt(i);
                    out << " " << orig_col << ":" << val;
                }
                out << "\n";
            }
        }

        out.close();
        std::cout << "Exported incidence matrix to " << filename << std::endl;
    }
};

#endif
