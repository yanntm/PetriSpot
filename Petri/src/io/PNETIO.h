/*
 * PNETIO.h
 *
 * PNET: binary container for a P/T net, the tool-to-tool net format
 * (INTEROP.md section 3). A 16-byte header followed by three KERS blocks:
 * flowPT (P x T), flowTP (P x T) and the initial marking (P x 1). Places
 * and transitions are identified by index; the loaded net names them
 * p<i> and t<i>.
 *
 * Header, little-endian:
 *   magic   : 4 bytes  "PNET"
 *   version : 1 byte   (= 1)
 *   flags   : 1 byte   (reserved, = 0)
 *   places  : 4 bytes  uint32
 *   trans   : 4 bytes  uint32
 *   padding : 2 bytes  (zero)
 */
#ifndef PETRI_IO_PNETIO_H_
#define PETRI_IO_PNETIO_H_

#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "core/MatrixCol.h"
#include "core/SparseArray.h"
#include "core/SparsePetriNet.h"
#include "io/SparseMatrixIO.h"

static const char PNET_MAGIC[4] = { 'P', 'N', 'E', 'T' };

template<typename T>
  class PNETIO
  {
    using IO = SparseMatrixIO<T>;

    static void expect (uint32_t got, uint32_t want, const std::string &what, const std::string &label)
    {
      if (got != want) {
        throw std::string ("PNET " + label + ": " + what + " has " + std::to_string (got) + " but the header says "
            + std::to_string (want));
      }
    }

  public:
    /** Write the net; false when the stream fails. */
    static bool write (const SparsePetriNet<T> &net, std::ostream &os)
    {
      uint32_t places = static_cast<uint32_t> (net.getPlaceCount ());
      uint32_t trans = static_cast<uint32_t> (net.getTransitionCount ());
      os.write (PNET_MAGIC, 4);
      os.put (1); // version
      os.put (0); // flags
      IO::writeLE (os, places);
      IO::writeLE (os, trans);
      os.put (0);
      os.put (0);
      if (!os) return false;
      if (!IO::write (net.getFlowPT (), os)) return false;
      if (!IO::write (net.getFlowTP (), os)) return false;
      MatrixCol<T> marking (places, 0);
      marking.appendColumn (SparseArray<T> (net.getMarks ()));
      return IO::write (marking, os);
    }

    /** Write a PNET file; errors on stderr, false on failure. */
    static bool write (const SparsePetriNet<T> &net, const std::string &filename)
    {
      std::ofstream ofs (filename, std::ios::binary | std::ios::out | std::ios::trunc);
      if (!ofs) {
        std::cerr << "Error: cannot open '" << filename << "' for writing\n";
        return false;
      }
      if (!write (net, ofs)) {
        std::cerr << "Error: write failed on '" << filename << "'\n";
        return false;
      }
      return true;
    }

    /** Read a net; throws std::string on any malformation. */
    static SparsePetriNet<T>* read (std::istream &is, const std::string &label)
    {
      uint8_t header[16];
      if (!is.read (reinterpret_cast<char*> (header), 16)) throw std::string ("Truncated PNET header in " + label);
      if (std::memcmp (header, PNET_MAGIC, 4) != 0) throw std::string ("Bad PNET magic in " + label);
      if (header[4] != 1) throw std::string ("Unsupported PNET version " + std::to_string (header[4]) + " in " + label);
      if (header[5] != 0) throw std::string ("Unsupported PNET flags in " + label);
      uint32_t places, trans;
      std::memcpy (&places, header + 6, 4);
      std::memcpy (&trans, header + 10, 4);
      places = IO::toLittleEndian (places);
      trans = IO::toLittleEndian (trans);

      uint32_t nr, nc;
      MatrixCol<T> pt = IO::read (is, nr, nc, label);
      expect (nr, places, "flowPT rows", label);
      expect (nc, trans, "flowPT columns", label);
      MatrixCol<T> tp = IO::read (is, nr, nc, label);
      expect (nr, places, "flowTP rows", label);
      expect (nc, trans, "flowTP columns", label);
      MatrixCol<T> m0 = IO::read (is, nr, nc, label);
      expect (nr, places, "marking rows", label);
      expect (nc, 1, "marking columns", label);
      std::vector<T> marks (places, 0);
      const SparseArray<T> &col = m0.getColumn (0);
      for (size_t i = 0; i < col.size (); ++i) marks[col.keyAt (i)] = col.valueAt (i);
      return new SparsePetriNet<T> (std::move (pt), std::move (tp), std::move (marks));
    }

    /** Read a PNET file; throws std::string on error. */
    static SparsePetriNet<T>* read (const std::string &filename)
    {
      std::ifstream ifs (filename, std::ios::binary);
      if (!ifs) throw std::string ("Cannot open PNET file: " + filename);
      return read (ifs, "'" + filename + "'");
    }
  };

#endif /* PETRI_IO_PNETIO_H_ */
