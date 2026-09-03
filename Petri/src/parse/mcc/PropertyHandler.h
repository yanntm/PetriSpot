/*
 * PropertyHandler.h
 *
 * expat SAX handler for the MCC property XML (property-set / property /
 * formula). Builds expr::Property values over a SparsePetriNet: places and
 * transitions are resolved by id, tokens-count becomes a sum of places,
 * is-fireable(t) is desugared into the conjunction of "m[p] >= w" over the
 * pre-arcs of t. Only the reachability fragment (EF phi, AG phi, EF deadlock)
 * yields a supported property; anything else is kept as Unsupported with a
 * comment.
 */
#ifndef PETRI_PARSE_MCC_PROPERTYHANDLER_H_
#define PETRI_PARSE_MCC_PROPERTYHANDLER_H_

#include <expat.h>
#include <string>
#include <unordered_map>
#include <vector>

#include "core/SparsePetriNet.h"
#include "expr/Property.h"

namespace petri::mcc
{

template<typename T>
  class PropertyHandler
  {
    using Expression = petri::expr::Expression;
    using LinearAtom = petri::expr::LinearAtom;
    using Property = petri::expr::Property;
    using PropertyKind = petri::expr::PropertyKind;

    /** An integer operand: a sum of places plus a constant. */
    struct IntOperand
    {
      std::vector<size_t> places;
      long long constant = 0;
    };

    /** One open XML element and the children completed so far. */
    struct Frame
    {
      std::string tag;
      std::string text;
      std::vector<Expression> bools;
      std::vector<IntOperand> ints;
      std::vector<size_t> places;
      std::vector<size_t> transitions;
      std::string temporal; // set by a temporal child: globally, finally...
      bool deadlock = false;
    };

    const SparsePetriNet<T> &net;
    std::unordered_map<std::string, size_t> placeIndex;
    std::unordered_map<std::string, size_t> transitionIndex;
    std::vector<Frame> stack;
    std::vector<Property> properties;
    Property current;
    bool failed = false; // current property left the supported fragment

  public:
    explicit PropertyHandler (const SparsePetriNet<T> &n)
        : net (n)
    {
      const auto &pn = net.getPnames ();
      for (size_t i = 0; i < pn.size (); ++i) placeIndex[pn[i]] = i;
      const auto &tn = net.getTnames ();
      for (size_t i = 0; i < tn.size (); ++i) transitionIndex[tn[i]] = i;
    }

    std::vector<Property>& getParseResult ()
    {
      return properties;
    }

    static void characters (void *userData, const XML_Char *chars, int length)
    {
      auto *self = static_cast<PropertyHandler<T>*> (userData);
      if (!self->stack.empty ()) {
        self->stack.back ().text.append (chars, static_cast<size_t> (length));
      }
    }

    static void startElement (void *userData, const XML_Char *name,
                              const XML_Char **)
    {
      auto *self = static_cast<PropertyHandler<T>*> (userData);
      std::string tag (name);
      if (tag == "property") {
        self->current = Property ();
        self->failed = false;
      }
      self->stack.push_back (Frame ());
      self->stack.back ().tag = std::move (tag);
    }

    static void endElement (void *userData, const XML_Char *name)
    {
      auto *self = static_cast<PropertyHandler<T>*> (userData);
      Frame f = std::move (self->stack.back ());
      self->stack.pop_back ();
      Frame *parent = self->stack.empty () ? nullptr : &self->stack.back ();
      self->close (f, parent, std::string (name));
    }

  private:
    static std::string trim (const std::string &s)
    {
      size_t b = s.find_first_not_of (" \t\r\n");
      if (b == std::string::npos) return "";
      size_t e = s.find_last_not_of (" \t\r\n");
      return s.substr (b, e - b + 1);
    }

    void unsupported (const std::string &why)
    {
      if (!failed) current.comment = why;
      failed = true;
    }

    /** Append the completed children of a closed frame to its parent. */
    static void forward (Frame &f, Frame *parent)
    {
      if (!parent) return;
      for (auto &b : f.bools) parent->bools.push_back (std::move (b));
      parent->deadlock = parent->deadlock || f.deadlock;
    }

    /** Build "sum(lhs) - sum(rhs) op rhs.constant - lhs.constant". */
    static Expression comparison (const IntOperand &lhs, const IntOperand &rhs,
                                  petri::expr::Cmp op)
    {
      LinearAtom a;
      for (size_t p : lhs.places) a.addTerm (p, 1);
      for (size_t p : rhs.places) a.addTerm (p, -1);
      a.constant = rhs.constant - lhs.constant;
      a.op = op;
      return Expression::makeAtom (std::move (a));
    }

    /** is-fireable(t): every pre-place holds at least the arc weight. */
    Expression fireable (size_t t) const
    {
      const auto &pre = net.getFlowPT ().getColumn (t);
      std::vector<Expression> kids;
      kids.reserve (pre.size ());
      for (size_t i = 0; i < pre.size (); ++i) {
        LinearAtom a;
        a.addTerm (pre.keyAt (i), 1);
        a.constant = static_cast<long long> (pre.valueAt (i));
        a.op = petri::expr::Cmp::GE;
        kids.push_back (Expression::makeAtom (std::move (a)));
      }
      if (kids.empty ()) return Expression::constant (true);
      if (kids.size () == 1) return std::move (kids[0]);
      return Expression::makeAnd (std::move (kids));
    }

    void pushBool (Frame *parent, Expression e)
    {
      if (parent) parent->bools.push_back (std::move (e));
    }

    void closeComparison (const Frame &f, Frame *parent, petri::expr::Cmp op)
    {
      if (f.ints.size () != 2) {
        unsupported ("comparison with " + std::to_string (f.ints.size ()) + " operands");
        pushBool (parent, Expression::constant (false));
        return;
      }
      pushBool (parent, comparison (f.ints[0], f.ints[1], op));
    }

    void close (Frame &f, Frame *parent, const std::string &tag)
    {
      using petri::expr::Cmp;
      if (tag == "property") {
        if (failed) current.kind = PropertyKind::Unsupported;
        properties.push_back (std::move (current));
      } else if (tag == "id") {
        current.name = trim (f.text);
      } else if (tag == "description" || tag == "property-set") {
        // ignored
      } else if (tag == "formula") {
        if (f.deadlock && current.kind != PropertyKind::Deadlock) {
          unsupported ("deadlock outside EF deadlock");
        } else if (f.bools.size () != 1 && current.kind != PropertyKind::Deadlock) {
          unsupported ("formula without a single body");
        }
        if (!f.bools.empty ()) current.body = std::move (f.bools[0]);
      } else if (tag == "all-paths" || tag == "exists-path") {
        bool isAll = tag == "all-paths";
        if (!parent || parent->tag != "formula") {
          unsupported ("nested path quantifier (CTL)");
        } else if (isAll && f.temporal == "globally") {
          current.kind = PropertyKind::Invariant;
        } else if (!isAll && f.temporal == "finally" && f.deadlock) {
          current.kind = PropertyKind::Deadlock;
        } else if (!isAll && f.temporal == "finally") {
          current.kind = PropertyKind::Reachability;
        } else {
          unsupported (tag + " " + f.temporal + " is not a reachability form");
        }
        if (current.kind == PropertyKind::Deadlock) f.deadlock = false;
        forward (f, parent);
      } else if (tag == "globally" || tag == "finally" || tag == "next"
          || tag == "until") {
        if (!parent || (parent->tag != "all-paths" && parent->tag != "exists-path")) {
          unsupported ("temporal operator " + tag + " outside a path quantifier (LTL)");
        }
        if (parent) parent->temporal = tag;
        forward (f, parent);
      } else if (tag == "deadlock") {
        if (parent) parent->deadlock = true;
      } else if (tag == "negation") {
        if (f.bools.size () == 1) pushBool (parent, Expression::makeNot (std::move (f.bools[0])));
        else unsupported ("negation arity");
        if (parent) parent->deadlock = parent->deadlock || f.deadlock;
      } else if (tag == "conjunction") {
        pushBool (parent, Expression::makeAnd (std::move (f.bools)));
      } else if (tag == "disjunction") {
        pushBool (parent, Expression::makeOr (std::move (f.bools)));
      } else if (tag == "integer-le") {
        closeComparison (f, parent, Cmp::LE);
      } else if (tag == "integer-lt") {
        closeComparison (f, parent, Cmp::LT);
      } else if (tag == "integer-ge") {
        closeComparison (f, parent, Cmp::GE);
      } else if (tag == "integer-gt") {
        closeComparison (f, parent, Cmp::GT);
      } else if (tag == "integer-eq") {
        closeComparison (f, parent, Cmp::EQ);
      } else if (tag == "integer-ne") {
        closeComparison (f, parent, Cmp::NEQ);
      } else if (tag == "tokens-count") {
        IntOperand op;
        op.places = std::move (f.places);
        if (parent) parent->ints.push_back (std::move (op));
      } else if (tag == "integer-constant") {
        IntOperand op;
        op.constant = std::stoll (trim (f.text));
        if (parent) parent->ints.push_back (std::move (op));
      } else if (tag == "place") {
        auto it = placeIndex.find (trim (f.text));
        if (it == placeIndex.end ()) {
          throw "Unknown place in property: " + trim (f.text);
        }
        if (parent) parent->places.push_back (it->second);
      } else if (tag == "transition") {
        auto it = transitionIndex.find (trim (f.text));
        if (it == transitionIndex.end ()) {
          throw "Unknown transition in property: " + trim (f.text);
        }
        if (parent) parent->transitions.push_back (it->second);
      } else if (tag == "is-fireable") {
        std::vector<Expression> kids;
        for (size_t t : f.transitions) kids.push_back (fireable (t));
        if (kids.empty ()) pushBool (parent, Expression::constant (false));
        else if (kids.size () == 1) pushBool (parent, std::move (kids[0]));
        else pushBool (parent, Expression::makeOr (std::move (kids)));
      } else if (tag == "place-bound") {
        unsupported ("place-bound (UpperBounds)");
        if (parent) parent->ints.push_back (IntOperand ());
      } else {
        unsupported ("unknown element " + tag);
      }
    }
  };

} // namespace petri::mcc

#endif /* PETRI_PARSE_MCC_PROPERTYHANDLER_H_ */
