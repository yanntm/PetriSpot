/*
 * PropertyHandler.h
 *
 * expat SAX handler for the MCC property XML (property-set / property /
 * formula). Builds expr::Property values over a SparsePetriNet: places and
 * transitions are resolved by id through NetResolver, tokens-count becomes a
 * sum of places, is-fireable is desugared into marking conditions. Every
 * boolean element yields a CTL formula node (a state predicate when it has
 * no temporal operator); the path quantifiers all-paths / exists-path with
 * globally, finally, next or until (before, reach) build the temporal nodes.
 * A closed formula is classified: EF phi / AG phi over a state predicate and
 * EF deadlock are the reachability kinds, place-bound (UpperBounds) a Bound,
 * anything else with a temporal operator a CTL property. Unknown elements
 * leave the property Unsupported with a comment.
 */
#ifndef PETRI_PARSE_MCC_PROPERTYHANDLER_H_
#define PETRI_PARSE_MCC_PROPERTYHANDLER_H_

#include <expat.h>
#include <string>
#include <vector>

#include "core/SparsePetriNet.h"
#include "expr/CtlFormula.h"
#include "expr/Property.h"
#include "parse/NetResolver.h"

namespace petri::mcc
{

template<typename T>
  class PropertyHandler
  {
    using Expression = petri::expr::Expression;
    using LinearAtom = petri::expr::LinearAtom;
    using Property = petri::expr::Property;
    using PropertyKind = petri::expr::PropertyKind;
    using CtlFormula = petri::expr::CtlFormula;
    using CtlOp = petri::expr::CtlOp;

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
      std::vector<CtlFormula> formulas; // boolean children, state or temporal
      std::vector<IntOperand> ints;
      std::vector<size_t> places;
      std::vector<size_t> transitions;
      std::string temporal; // set by a temporal child: globally, finally, next, until
      bool bound = false;   // a place-bound child (UpperBounds)
    };

    NetResolver<T> net;
    std::vector<Frame> stack;
    std::vector<Property> properties;
    Property current;
    bool failed = false; // current property left the supported fragment

  public:
    explicit PropertyHandler (const SparsePetriNet<T> &n)
        : net (n)
    {
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

    /** Append the completed formulas of a closed frame to its parent. */
    static void forward (Frame &f, Frame *parent)
    {
      if (!parent) return;
      for (auto &b : f.formulas) parent->formulas.push_back (std::move (b));
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
      a.normalize ();
      return Expression::makeAtom (std::move (a));
    }

    void push (Frame *parent, CtlFormula f)
    {
      if (parent) parent->formulas.push_back (std::move (f));
    }

    void pushPred (Frame *parent, Expression e)
    {
      push (parent, CtlFormula::predicate (std::move (e)));
    }

    void closeComparison (const Frame &f, Frame *parent, petri::expr::Cmp op)
    {
      if (f.ints.size () != 2) {
        unsupported ("comparison with " + std::to_string (f.ints.size ()) + " operands");
        pushPred (parent, Expression::constant (false));
        return;
      }
      pushPred (parent, comparison (f.ints[0], f.ints[1], op));
    }

    static bool allPreds (const std::vector<CtlFormula> &fs)
    {
      for (const auto &f : fs) if (f.op != CtlOp::Pred) return false;
      return true;
    }

    /** and / or over children: one predicate when every child is one, a CTL node otherwise. */
    void closeBool (Frame &f, Frame *parent, CtlOp op)
    {
      if (allPreds (f.formulas)) {
        std::vector<Expression> es;
        for (auto &k : f.formulas) es.push_back (std::move (k.pred));
        pushPred (parent, op == CtlOp::And ? Expression::makeAnd (std::move (es)) : Expression::makeOr (std::move (es)));
      } else {
        push (parent, CtlFormula::nary (op, std::move (f.formulas)));
      }
    }

    /** all-paths / exists-path: the temporal child names the operator. */
    void closeQuantifier (Frame &f, Frame *parent, bool isAll)
    {
      const std::string &t = f.temporal;
      CtlOp op;
      size_t arity = 1;
      if (t == "globally") op = isAll ? CtlOp::AG : CtlOp::EG;
      else if (t == "finally") op = isAll ? CtlOp::AF : CtlOp::EF;
      else if (t == "next") op = isAll ? CtlOp::AX : CtlOp::EX;
      else if (t == "until") { op = isAll ? CtlOp::AU : CtlOp::EU; arity = 2; }
      else {
        unsupported ("path quantifier without a temporal operator");
        return;
      }
      if (f.formulas.size () != arity) {
        unsupported (t + " with " + std::to_string (f.formulas.size ()) + " operands");
        return;
      }
      if (arity == 1) push (parent, CtlFormula::unary (op, std::move (f.formulas[0])));
      else push (parent, CtlFormula::binary (op, std::move (f.formulas[0]), std::move (f.formulas[1])));
    }

    /** The closed formula: reachability fragment, or CTL. */
    void classify (CtlFormula f)
    {
      if (f.op == CtlOp::EF && f.kids[0].op == CtlOp::Pred) {
        current.kind = PropertyKind::Reachability;
        current.body = std::move (f.kids[0].pred);
      } else if (f.op == CtlOp::EF && f.kids[0].op == CtlOp::Deadlock) {
        current.kind = PropertyKind::Deadlock;
      } else if (f.op == CtlOp::AG && f.kids[0].op == CtlOp::Pred) {
        current.kind = PropertyKind::Invariant;
        current.body = std::move (f.kids[0].pred);
      } else {
        current.kind = PropertyKind::CTL;
        current.ctl = std::move (f);
      }
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
        if (f.bound) {
          if (f.ints.size () != 1 || !f.formulas.empty ()) {
            unsupported ("place-bound combined with other operators");
          } else {
            LinearAtom a;
            for (size_t p : f.ints[0].places) a.addTerm (p, 1);
            a.normalize ();
            a.op = Cmp::GE;
            current.kind = PropertyKind::Bound;
            current.body = Expression::makeAtom (std::move (a));
          }
        } else if (f.formulas.size () != 1) {
          unsupported ("formula without a single body");
        } else {
          classify (std::move (f.formulas[0]));
        }
      } else if (tag == "all-paths" || tag == "exists-path") {
        closeQuantifier (f, parent, tag == "all-paths");
      } else if (tag == "globally" || tag == "finally" || tag == "next" || tag == "until") {
        if (!parent || (parent->tag != "all-paths" && parent->tag != "exists-path")) {
          unsupported ("temporal operator " + tag + " outside a path quantifier (LTL)");
        }
        if (parent) parent->temporal = tag;
        forward (f, parent);
      } else if (tag == "before" || tag == "reach") {
        if (!parent || parent->tag != "until") unsupported (tag + " outside until");
        if (f.formulas.size () != 1) unsupported (tag + " without a single operand");
        forward (f, parent);
      } else if (tag == "deadlock") {
        push (parent, CtlFormula::leaf (CtlOp::Deadlock));
      } else if (tag == "negation") {
        if (f.formulas.size () != 1) {
          unsupported ("negation arity");
        } else if (f.formulas[0].op == CtlOp::Pred) {
          pushPred (parent, Expression::makeNot (std::move (f.formulas[0].pred)));
        } else {
          push (parent, CtlFormula::unary (CtlOp::Not, std::move (f.formulas[0])));
        }
      } else if (tag == "conjunction") {
        closeBool (f, parent, CtlOp::And);
      } else if (tag == "disjunction") {
        closeBool (f, parent, CtlOp::Or);
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
        auto p = net.place (trim (f.text));
        if (!p) throw "Unknown place in property: " + trim (f.text);
        if (parent) parent->places.push_back (*p);
      } else if (tag == "transition") {
        auto t = net.transition (trim (f.text));
        if (!t) throw "Unknown transition in property: " + trim (f.text);
        if (parent) parent->transitions.push_back (*t);
      } else if (tag == "is-fireable") {
        pushPred (parent, net.anyFireable (f.transitions));
      } else if (tag == "place-bound") {
        IntOperand op;
        op.places = std::move (f.places);
        if (parent) {
          parent->ints.push_back (std::move (op));
          parent->bound = true;
        }
      } else {
        unsupported ("unknown element " + tag);
      }
    }
  };

} // namespace petri::mcc

#endif /* PETRI_PARSE_MCC_PROPERTYHANDLER_H_ */
