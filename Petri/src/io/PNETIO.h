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
 *
 * After the three mandatory blocks a net may carry optional *named* blocks,
 * each one a 12-byte framing followed by an ordinary KERS payload:
 *
 *   name    : 8 bytes  ASCII, zero-padded (e.g. "TMULT")
 *   length  : 4 bytes  uint32, the payload's byte count
 *
 * They are read only by a caller that asks for them, an unknown name is
 * skipped by its length, and a producer emits only what it has: absence of a
 * block is how "not available" is said. Extensibility is by new names, so the
 * framing never changes; see HSC_PLAN.md sections 10 to 13 for the blocks in
 * use (`TMULT`, transition multiplicities).
 */
#ifndef PETRI_IO_PNETIO_H_
#define PETRI_IO_PNETIO_H_

#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
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

  public:
    /** The optional named blocks, in file order. */
    using Blocks = std::vector<std::pair<std::string, MatrixCol<T>>>;

    /** The block of that name, or nullptr. */
    static const MatrixCol<T>* find (const Blocks &blocks, const std::string &name)
    {
      for (const auto &b : blocks)
        if (b.first == name) return &b.second;
      return nullptr;
    }

  private:
    /** A block name as 8 zero-padded bytes; throws when it does not fit. */
    static void writeName (std::ostream &os, const std::string &name)
    {
      if (name.empty () || name.size () > 8) {
        throw std::string ("PNET block name '" + name + "' must be 1 to 8 characters");
      }
      char buf[8] = { 0 };
      std::memcpy (buf, name.data (), name.size ());
      os.write (buf, 8);
    }

    static std::string trimName (const char *raw)
    {
      std::string name (raw, 8);
      while (!name.empty () && (name.back () == '\0' || name.back () == ' ')) name.pop_back ();
      return name;
    }

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

    /** Write the net followed by \p blocks, in the order given. */
    static bool write (const SparsePetriNet<T> &net, std::ostream &os, const Blocks &blocks)
    {
      if (!write (net, os)) return false;
      for (const auto &b : blocks) {
        std::ostringstream payload (std::ios::binary);
        if (!IO::write (b.second, payload)) return false;
        const std::string bytes = payload.str ();
        writeName (os, b.first);
        IO::writeLE (os, static_cast<uint32_t> (bytes.size ()));
        os.write (bytes.data (), static_cast<std::streamsize> (bytes.size ()));
        if (!os) return false;
      }
      return true;
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

    /** Read a net; throws std::string on any malformation. \p blocks, when
     * given, receives the optional named blocks that follow (see the header
     * comment); without it they are simply not read. */
    static SparsePetriNet<T>* read (std::istream &is, const std::string &label,
                                    Blocks *blocks = nullptr)
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
      if (blocks != nullptr) readBlocks (is, label, *blocks);
      return new SparsePetriNet<T> (std::move (pt), std::move (tp), std::move (marks));
    }

    /** Read the optional named blocks to the end of the stream. */
    static void readBlocks (std::istream &is, const std::string &label, Blocks &blocks)
    {
      for (;;) {
        char raw[8];
        is.read (raw, 8);
        if (is.gcount () == 0) return;  // clean end of file
        if (is.gcount () != 8) throw std::string ("Truncated PNET block name in " + label);
        uint32_t len = 0;
        if (!is.read (reinterpret_cast<char*> (&len), 4)) {
          throw std::string ("Truncated PNET block length in " + label);
        }
        len = IO::toLittleEndian (len);
        std::string payload (len, '\0');
        if (len > 0 && !is.read (payload.data (), len)) {
          throw std::string ("Truncated PNET block payload in " + label);
        }
        std::istringstream ss (payload, std::ios::binary);
        uint32_t nr, nc;
        MatrixCol<T> m = IO::read (ss, nr, nc, label);
        blocks.emplace_back (trimName (raw), std::move (m));
      }
    }

    /** Read a PNET file; throws std::string on error. */
    static SparsePetriNet<T>* read (const std::string &filename, Blocks *blocks = nullptr)
    {
      std::ifstream ifs (filename, std::ios::binary);
      if (!ifs) throw std::string ("Cannot open PNET file: " + filename);
      return read (ifs, "'" + filename + "'", blocks);
    }
  };

#endif /* PETRI_IO_PNETIO_H_ */
