// FIXME Licence GPLv3, ask @ytm?

#pragma once

#include <sstream>
#include <vector>
#include "SparseIntArray.h"

#include <spot/kripke/kripke.hh>


struct SparseIntArray_equal
{
  bool operator()(const SparseIntArray & lhs, const SparseIntArray & rhs) const
  {
    return lhs == rhs;
  }
};

struct SparseIntArray_hash
{
  size_t operator()(const SparseIntArray & that) const
  {
    return that.hash();
  }
};


class PT_iterator
{
 public:

 PT_iterator(const SparseIntArray & src, const MatrixCol & combFlow, int * list, spot::cube cond)
   : src_(src),combFlow_(combFlow),list_(list), i_(0), cond_(cond)
    { }

  ~PT_iterator()
    {
      // FIXME can we delete elements in succs_ since they
      // will be inserted inside of the hashmap?
      // Since Spot relies on automatic memory management
      // we must add methods to delete states
	  delete[] list_;
    }

  void next()
  {
    ++i_;
  }

  bool done()
  {
    return i_ == list_[0];
  }

  SparseIntArray state()
  {
	  //dropEmpty(list);
	  return  SparseIntArray::sumProd(1, src_, 1, combFlow_.getColumn(list_[i_+1]));
  }

  spot::cube condition()
    {
      return cond_;
    }
 private:
  SparseIntArray src_;
  const MatrixCol & combFlow_;
  int* list_;
  unsigned i_;
  spot::cube cond_;
};

template<>
class spot::kripkecube<SparseIntArray, PT_iterator> final
{
 public:

  kripkecube(const SparsePetriNet& sr, const std::vector<std::string> & observed_aps, petri::expr::AtomicPropManager & apm)
    : sr_(sr), observed_aps_(observed_aps), combFlow(sr.getPnames().size(), 0)
  {
    // FIXME : do we really need this?
    for (int i = 0 ;  i < sr.getFlowPT().getColumnCount() ; i ++)
      {
	combFlow.appendColumn(SparseIntArray::sumProd(-1, sr.getFlowPT().getColumn(i),
						      1, sr.getFlowTP().getColumn(i)));
      }

    for (auto & ap : observed_aps) {
    	int index = atoi( ap.c_str() + 1);
    	aps_.push_back(apm.getAtoms().at(index).second);
    }

    // Check that the walker fits the requirements
    static_assert(spot::is_a_kripkecube_ptr<decltype(this),
		  SparseIntArray, PT_iterator>::value,
		  "error: does not match the kripkecube requirements");
  }

  ~kripkecube() = default;


  // ------------------------------------------------------------------------
  // Match Spot Kripkecube "interface"
  // ------------------------------------------------------------------------

  /// \brief Returns the initial SparseIntArray of the System. The \a tid parameter
  /// is used internally for sharing this structure among threads.
  SparseIntArray initial(unsigned tid = 0)
  {
    (void) tid;

    SparseIntArray state(sr_.getMarks());
    return state;
  }

  /// \brief Returns the number of threads that are handled by the kripkecube
  unsigned get_threads()
  {
    return 8; // FIXME
  }

  /// \brief Provides a string representation of the parameter state
  std::string to_string(const SparseIntArray & s, unsigned tid = 0) const
    {
      (void) tid;

      std::ostringstream stream;
      s.print(stream);
      return stream.str();
    }

  /// \brief Returns an iterator over the successors of the parameter state.
  PT_iterator* succ(const SparseIntArray & s, unsigned tid = 0)
  {
    (void) tid;
    int* list = computeEnabled(s);

    auto cs = spot::cubeset(observed_aps_.size());
    cube cond = cs.alloc();
    int i =0;
    for (auto & ap : aps_) {
    	if (ap->eval(s)) {
    		cs.set_true_var(cond, i);
    	} else {
    		cs.set_false_var(cond, i);
    	}
    	++i;
    }

    return new PT_iterator(s,combFlow,list,cond);
  }

  /// \brief Allocation and deallocation of iterator is costly. This
  /// method allows to reuse old iterators.
  void recycle(PT_iterator* it, unsigned tid)
  {
    (void) tid;

    delete it; // FIXME use tid to build a cache per thread
  }

  /// \brief This method returns the observed atomic propositions
  const std::vector<std::string> ap()
  {
    return observed_aps_;
  }

  // ------------------------------------------------------------------------
  // End matching the  Spot Kripkecube "interface"
  // ------------------------------------------------------------------------

 private:
  int* computeEnabled(const SparseIntArray& state)
  {
    int * list  = new int [sr_.getTnames().size()+1];
    memset(list,0, (sr_.getTnames().size()+1) * sizeof(int));

    int li = 1;
    int t = 0;
    for (int e = sr_.getTnames().size(); t < e; t++)
      {
	if (SparseIntArray::greaterOrEqual(state, sr_.getFlowPT().getColumn(t)))
	  {
	    list[li++] = t;
	  }
      }
    list[0] = li -1 ;
    return list;
  }

  void dropEmpty(int* enabled)
  {
    for (int i = enabled[0] ; i  >= 1  ; i--)
      {
	int t = enabled [i];
	if (combFlow.getColumn(t).size() == 0)
	  {
	    dropAt(enabled,i);
	  }
      }
  }

  void dropAt (int* enabled, int index)
  {
    if (index < enabled[0])
      {
	enabled[index] = enabled[enabled[0]];
      }
    enabled [0] --;
  }

  SparseIntArray fire (int t, const SparseIntArray & state)
  {
    // NB no enabling check
    return SparseIntArray::sumProd(1, state, 1, combFlow.getColumn(t));
  }

 private:
  const SparsePetriNet& sr_;
  std::vector<std::string> observed_aps_;
  std::vector<petri::expr::Expression *> aps_;
  MatrixCol combFlow; // FIXME one per thread?
};

typedef spot::kripkecube<SparseIntArray, PT_iterator> Petricube;
typedef spot::kripkecube<SparseIntArray, PT_iterator>* Petricube_ptr;
