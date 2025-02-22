#ifndef PTNETHANDLER_H_
#define PTNETHANDLER_H_

#include "SparsePetriNet.h"
#include "ddd/util/ext_hash_map.hh"
#include <expat.h>
#include <stack>
#include <iostream>

template<typename T>
  class PTNetHandler
  {
    // Enable for verbose debugging output during development.
    static inline const bool DEBUG = false;

    std::stack<void*> stack;
    SparsePetriNet<T> *net;
    typedef std::pair<bool, int> node_t;
    typedef ext_hash_map<std::string, node_t> index_t;
    index_t index;
    typedef std::pair<std::pair<std::string, std::string>, T> arc_t;
    typedef std::vector<arc_t*> arcs_t;
    arcs_t topatch;

    std::string lastseen;
    bool readtext;
    long lastint;
    bool readint;
    bool inOpaqueToolSpecific;
    bool doIt;
    std::string textBuffer; // New buffer to accumulate character data

  public:
    PTNetHandler ()
        : net (new SparsePetriNet<T> ()), readtext (false), lastint (-1), readint (
            false), inOpaqueToolSpecific (false), doIt (false), textBuffer ("")
    {
      if (DEBUG) std::cout << "PTNetHandler initialized" << std::endl;
    }

    static void characters (void *userData, const XML_Char *chars, int length)
    {
      PTNetHandler<T> *tthis = (PTNetHandler<T>*) userData;
      if (tthis->inOpaqueToolSpecific) {
        return;
      }
      std::string data (chars, length);
      if (tthis->readtext || tthis->readint) { // Accumulate regardless of doIt for now
        tthis->textBuffer += data; // Append to buffer
        if (DEBUG) std::cout << "Accumulating characters: '" << data
            << "' -> buffer = '" << tthis->textBuffer << "'" << std::endl;
      } else if (DEBUG) {
        std::cout << "Characters ignored: '" << data << "'" << std::endl;
      }
    }

    static void startElement (void *userData, const XML_Char *name,
                              const XML_Char **atts)
    {
      PTNetHandler<T> *tthis = (PTNetHandler<T>*) userData;
      if (tthis->inOpaqueToolSpecific) {
        return;
      }

      std::string baliseName (name);
      if (DEBUG) std::cout << "Start element: <" << baliseName << ">"
          << std::endl;

      // Clear textBuffer at the start of a new element we care about
      if ("name" == baliseName || "initialMarking" == baliseName
          || "inscription" == baliseName) {
        tthis->textBuffer.clear ();
        if (DEBUG) std::cout << "  Cleared textBuffer for " << baliseName
            << std::endl;
      }

      if ("net" == baliseName) {
        for (int i = 0; atts[i] != nullptr; i += 2) {
          std::string bname = atts[i];
          std::string bval = atts[i + 1];
          if (DEBUG) std::cout << "  Attribute: " << bname << " = " << bval
              << std::endl;
          if ("id" == bname) {
            tthis->net->setName (bval);
          } else if ("type" == bname) {
            if ("http://www.pnml.org/version-2009/grammar/ptnet" != bval) {
              throw "Net is not a P/T net-> Colors are not supported currently.";
            }
          }
        }
        tthis->stack.push (&tthis->net);
      } else if ("name" == baliseName) {
        tthis->readtext = true;
        if (DEBUG) std::cout << "  Expecting text for name" << std::endl;
      } else if ("place" == baliseName) {
        SparsePetriNet<T> *pn = tthis->net;
        std::string id;
        for (int i = 0; atts[i] != nullptr; i += 2) {
          std::string bname = atts[i];
          std::string bval = atts[i + 1];
          if ("id" == bname) {
            id = bval;
            break;
          }
        }
        size_t pid = pn->addPlace (id, 0);
        index_t::accessor acc;
        tthis->index.insert (acc, id);
        acc->second =
          { true, static_cast<int> (pid) };
        tthis->stack.push (reinterpret_cast<void*> (pid));
        if (DEBUG) std::cout << "  Added place '" << id << "' with pid " << pid
            << std::endl;
      } else if ("initialMarking" == baliseName) {
        tthis->readint = true;
        if (DEBUG) std::cout << "  Expecting int for initialMarking"
            << std::endl;
      } else if ("inscription" == baliseName) {
        tthis->readint = true;
        if (DEBUG) std::cout << "  Expecting int for inscription" << std::endl;
      } else if ("transition" == baliseName) {
        SparsePetriNet<T> *pn = tthis->net;
        std::string id;
        for (int i = 0; atts[i] != nullptr; i += 2) {
          std::string bname = atts[i];
          std::string bval = atts[i + 1];
          if ("id" == bname) {
            id = bval;
            break;
          }
        }
        size_t tid = pn->addTransition (id);
        index_t::accessor acc;
        tthis->index.insert (acc, id);
        acc->second =
          { false, static_cast<int> (tid) };
        tthis->stack.push (reinterpret_cast<void*> (tid));
        if (DEBUG) std::cout << "  Added transition '" << id << "' with tid "
            << tid << std::endl;
      } else if ("arc" == baliseName) {
        std::string source, target;
        for (int i = 0; atts[i] != nullptr; i += 2) {
          std::string bname = atts[i];
          std::string bval = atts[i + 1];
          if ("source" == bname) source = bval;
          else if ("target" == bname) target = bval;
        }
        arc_t *arc = new arc_t (
          {
            { source, target }, 1 });
        tthis->stack.push (arc);
        if (DEBUG) std::cout << "  Added arc " << source << " -> " << target
            << " (weight=1)" << std::endl;
      } else if ("toolspecific" == baliseName) {
        tthis->inOpaqueToolSpecific = true;
        if (DEBUG) std::cout << "  Entering toolspecific section" << std::endl;
      } else if ("text" == baliseName) {
        tthis->doIt = true;
        if (DEBUG) std::cout << "  Enabling text/int parsing" << std::endl;
      } else if ("graphics" == baliseName || "offset" == baliseName
          || "position" == baliseName || "fill" == baliseName
          || "line" == baliseName || "dimension" == baliseName) {
        // Skip
      } else if ("pnml" == baliseName) {
        // Skip
      } else if ("page" == baliseName) {
        // Skip
      } else {
        std::cerr << "Unknown XML tag in source file: " << baliseName
            << std::endl;
      }
    }

    static void endElement (void *userData, const XML_Char *name)
    {
      PTNetHandler<T> *tthis = (PTNetHandler<T>*) userData;
      std::string baliseName (name);
      if (DEBUG) std::cout << "End element: </" << baliseName << ">"
          << std::endl;

      if ("toolspecific" == baliseName) {
        tthis->inOpaqueToolSpecific = false;
        if (DEBUG) std::cout << "  Exiting toolspecific section" << std::endl;
      } else if (tthis->inOpaqueToolSpecific) {
        return;
      } else if ("net" == baliseName) {
        tthis->stack.pop ();
        assert(tthis->stack.empty());
      } else if ("name" == baliseName) {
        tthis->lastseen = tthis->textBuffer; // Assign accumulated text
        tthis->readtext = false;
        tthis->textBuffer.clear ();
        if (DEBUG) std::cout << "  Set lastseen = '" << tthis->lastseen
            << "' and cleared textBuffer" << std::endl;
      } else if ("place" == baliseName) {
        tthis->stack.pop ();
      } else if ("transition" == baliseName) {
        tthis->stack.pop ();
      } else if ("arc" == baliseName) {
        arc_t *arc = static_cast<arc_t*> (tthis->stack.top ());
        index_t::const_accessor accsrc, acctgt;
        bool oksrc = tthis->index.find (accsrc, arc->first.first);
        bool oktgt = tthis->index.find (acctgt, arc->first.second);
        if (oktgt && oksrc) {
          if (accsrc->second.first) {
            tthis->net->addPreArc (accsrc->second.second, acctgt->second.second,
                                   arc->second);
            if (DEBUG) std::cout << "  Added pre-arc P" << accsrc->second.second
                << " -> T" << acctgt->second.second << " (weight="
                << arc->second << ")" << std::endl;
          } else {
            tthis->net->addPostArc (acctgt->second.second,
                                    accsrc->second.second, arc->second);
            if (DEBUG) std::cout << "  Added post-arc T"
                << accsrc->second.second << " -> P" << acctgt->second.second
                << " (weight=" << arc->second << ")" << std::endl;
          }
          delete arc;
        } else {
          tthis->topatch.push_back (arc);
          if (DEBUG) std::cout << "  Deferred arc " << arc->first.first
              << " -> " << arc->first.second << std::endl;
        }
        tthis->stack.pop ();
      } else if ("text" == baliseName) {
        tthis->doIt = false;
        if (DEBUG) std::cout << "  Disabling text/int parsing" << std::endl;
      } else if ("initialMarking" == baliseName) {
        size_t p = reinterpret_cast<size_t> (tthis->stack.top ());
        try {
          tthis->lastint = std::stol (tthis->textBuffer); // Convert accumulated text to int
          if (DEBUG) std::cout << "  Parsed initialMarking from '"
              << tthis->textBuffer << "' -> lastint = " << tthis->lastint
              << std::endl;
        } catch (const std::exception &e) {
          if (DEBUG) std::cout << "  Failed to parse initialMarking from '"
              << tthis->textBuffer << "': " << e.what () << std::endl;
          tthis->lastint = 0; // Default to 0 on error
        }
        tthis->net->setMarking (p, tthis->lastint);
        if (DEBUG) std::cout << "Initial marking for place " << p << " = "
            << tthis->net->getPnames ()[p] << " is " << tthis->lastint
            << std::endl;
        tthis->readint = false;
        tthis->textBuffer.clear ();
        tthis->lastint = -1;
        if (DEBUG) std::cout
            << "  Reset readint, cleared textBuffer, and reset lastint"
            << std::endl;
      } else if ("inscription" == baliseName) {
        arc_t *arc = static_cast<arc_t*> (tthis->stack.top ());
        try {
          tthis->lastint = std::stol (tthis->textBuffer); // Convert accumulated text to int
          if (DEBUG) std::cout << "  Parsed inscription from '"
              << tthis->textBuffer << "' -> lastint = " << tthis->lastint
              << std::endl;
        } catch (const std::exception &e) {
          if (DEBUG) std::cout << "  Failed to parse inscription from '"
              << tthis->textBuffer << "': " << e.what () << std::endl;
          tthis->lastint = 1; // Default to 1 on error (arc weight)
        }
        arc->second = tthis->lastint;
        tthis->readint = false;
        tthis->textBuffer.clear ();
        tthis->lastint = -1;
        if (DEBUG) std::cout << "  Set arc weight to " << arc->second
            << ", reset readint, cleared textBuffer, and reset lastint"
            << std::endl;
      } else if ("pnml" == baliseName) {
        for (arc_t *&arc : tthis->topatch) {
          index_t::const_accessor accsrc, acctgt;
          bool oksrc = tthis->index.find (accsrc, arc->first.first);
          bool oktgt = tthis->index.find (acctgt, arc->first.second);
          if (oktgt && oksrc) {
            if (accsrc->second.first) {
              tthis->net->addPreArc (accsrc->second.second,
                                     acctgt->second.second, arc->second);
              if (DEBUG) std::cout << "  Patched pre-arc P"
                  << accsrc->second.second << " -> T" << acctgt->second.second
                  << " (weight=" << arc->second << ")" << std::endl;
            } else {
              tthis->net->addPostArc (acctgt->second.second,
                                      accsrc->second.second, arc->second);
              if (DEBUG) std::cout << "  Patched post-arc T"
                  << accsrc->second.second << " -> P" << acctgt->second.second
                  << " (weight=" << arc->second << ")" << std::endl;
            }
            delete arc;
          } else {
            std::string err =
                "Problem when linking arc : source or target node not found <"
                    + accsrc->first + "," + acctgt->first + ">";
            delete arc;
            throw err.c_str ();
          }
        }
        tthis->topatch.clear ();
        if (DEBUG) std::cout << "  Cleared topatch" << std::endl;
      } else if ("graphics" == baliseName || "offset" == baliseName
          || "position" == baliseName || "fill" == baliseName
          || "line" == baliseName || "dimension" == baliseName) {
        // Skip
      } else {
        std::cerr << "Unknown XML tag in source file: " << baliseName
            << std::endl;
      }
    }

    SparsePetriNet<T>* getParseResult ()
    {
      return net;
    }
  };

#endif
