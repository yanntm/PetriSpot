// kersconv: convert between KERS binary and ASCII sparse matrix formats.
//
// Usage:
//   kersconv --decode <input.kers> [output.txt]   (KERS → ASCII, default stdout)
//   kersconv --encode <input.txt>  <output.kers>  (ASCII → KERS)
//
// ASCII format (same as --exportAsMatrix):
//   <nrows> <ncols>
//   <rowIdx> <colIdx>:<val> ...     (one non-empty row per line)

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <numeric>
#include "MatrixCol.h"
#include "SparseArray.h"
#include "SparseMatrixIO.h"

// Decode KERS → ASCII
static int decode (const std::string &kersFile, std::ostream &out)
{
  MatrixCol<long> M = SparseMatrixIO<long>::read (kersFile);
  if (M.getColumnCount () == 0 && M.getRowCount () == 0)
    return 1;

  size_t rows = M.getRowCount ();
  size_t cols = M.getColumnCount ();
  out << rows << " " << cols << "\n";

  // Transpose so we can iterate row-by-row
  MatrixCol<long> T = M.transpose ();
  const auto &tcols = T.getColumns ();
  for (size_t r = 0; r < T.getColumnCount (); ++r) {
    const SparseArray<long> &row = tcols[r];
    if (row.size () == 0) continue;
    out << r;
    for (size_t k = 0; k < row.size (); ++k)
      out << " " << row.keyAt (k) << ":" << row.valueAt (k);
    out << "\n";
  }
  return 0;
}

// Encode ASCII → KERS
static int encode (const std::string &asciiFile, const std::string &kersFile)
{
  std::ifstream in (asciiFile);
  if (!in) {
    std::cerr << "Error: cannot open " << asciiFile << std::endl;
    return 1;
  }

  size_t nrows, ncols;
  in >> nrows >> ncols;
  if (!in) {
    std::cerr << "Error: bad header in " << asciiFile << std::endl;
    return 1;
  }

  MatrixCol<long> M (nrows, ncols);

  std::string line;
  std::getline (in, line); // consume rest of header line
  while (std::getline (in, line)) {
    if (line.empty ()) continue;
    std::istringstream ss (line);
    size_t rowIdx;
    ss >> rowIdx;
    std::string token;
    // We accumulate (col, val) per row, then set them via transpose trick.
    // Since MatrixCol is column-sparse, use a transposed builder.
    while (ss >> token) {
      auto sep = token.find (':');
      if (sep == std::string::npos) {
        std::cerr << "Error: malformed token '" << token << "'" << std::endl;
        return 1;
      }
      size_t colIdx = std::stoul (token.substr (0, sep));
      long val = std::stol (token.substr (sep + 1));
      // Column colIdx, row rowIdx
      M.getColumn (colIdx).append ((int)rowIdx, val);
    }
  }

  if (!SparseMatrixIO<long>::write (M, kersFile))
    return 1;
  std::cerr << "Encoded " << nrows << "x" << ncols << " matrix to " << kersFile << std::endl;
  return 0;
}

int main (int argc, char *argv[])
{
  if (argc < 3) {
    std::cerr << "Usage:\n"
              << "  kersconv --decode <input.kers> [output.txt]\n"
              << "  kersconv --encode <input.txt>  <output.kers>\n";
    return 1;
  }

  std::string mode = argv[1];

  if (mode == "--decode") {
    std::string kersFile = argv[2];
    if (argc >= 4) {
      std::ofstream out (argv[3]);
      if (!out) { std::cerr << "Error: cannot open " << argv[3] << std::endl; return 1; }
      return decode (kersFile, out);
    } else {
      return decode (kersFile, std::cout);
    }
  } else if (mode == "--encode") {
    if (argc < 4) { std::cerr << "Error: --encode requires input and output paths.\n"; return 1; }
    return encode (argv[2], argv[3]);
  } else {
    std::cerr << "Unknown mode: " << mode << std::endl;
    return 1;
  }
}
