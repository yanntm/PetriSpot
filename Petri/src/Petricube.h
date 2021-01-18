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

 PT_iterator(SparseIntArray* succs, unsigned size)
   : succs_(succs), i_(0), size_(size)
    { }

  ~PT_iterator()
    {
      // FIXME can we delete elements in succs_ since they
      // will be inserted inside of the hashmap?
      // Since Spot relies on automatic memory management
      // we must add methods to delete states
    }

  void next()
  {
    ++i_;
  }

  bool done()
  {
    return i_ == size_;
  }

  SparseIntArray state()
  {
    return succs_[i_];
  }

  spot::cube condition()
    {
      assert(false);
      return nullptr;
    }
 private:
  SparseIntArray* succs_;
  unsigned i_;
  unsigned size_;
};

template<>
class spot::kripkecube<SparseIntArray, PT_iterator> final
{
 public:

  kripkecube(const SparsePetriNet& sr, const std::vector<std::string> & observed_aps)
    : sr_(sr), observed_aps_(observed_aps), combFlow(sr.getPnames().size(), 0)
  {
    // FIXME : do we really need this?
    for (int i = 0 ;  i < sr.getFlowPT().getColumnCount() ; i ++)
      {
	combFlow.appendColumn(SparseIntArray::sumProd(-1, sr.getFlowPT().getColumn(i),
						      1, sr.getFlowTP().getColumn(i)));
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
    return 1; // FIXME
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
    dropEmpty(list);
    SparseIntArray* succ = new SparseIntArray[list[0]];
    for (int ti = 1 ; ti -1 < list[0] ; ti++)
      {
	succ[ti-1] = fire(list[ti], s /*FIXME is s the good parameter here?*/);
	//i++; FIXME useless here?
      }
    // FIXME can we delete here list?

    return new PT_iterator(succ, list[0]);
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

  MatrixCol combFlow; // FIXME one per thread?
};

typedef spot::kripkecube<SparseIntArray, PT_iterator> Petricube;
typedef spot::kripkecube<SparseIntArray, PT_iterator>* Petricube_ptr;
