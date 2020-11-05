// FIXME Licence GPLv3, ask @ytm?

#pragma once


#include <spot/kripke/kripke.hh>

class Petricube
{
 public:
  Petricube()
    {
    // Check that the walker fits the requirements
    static_assert(spot::is_a_kripkecube_ptr<decltype(this),
		  State, SuccIterator>::value,
		  "error: does not match the kripkecube requirements");
    }

  ~Petricube() = default;


  // ------------------------------------------------------------------------
  // Match Spot Kripkecube "interface"
  // ------------------------------------------------------------------------
  class PT_iterator; // Forward declaration
  using State = int; // FIXME
  using SuccIterator = PT_iterator;

  /// \brief Returns the initial State of the System. The \a tid parameter
  /// is used internally for sharing this structure among threads.
  State initial(unsigned tid)
  {
    (void) tid;
    assert(false);
    return 0;
  }

  /// \brief Returns the number of threads that are handled by the kripkecube
  unsigned get_threads()
  {
    assert(false);
    return 0;
  }

  /// \brief Provides a string representation of the parameter state
  std::string to_string(const State s, unsigned tid) const
    {
      (void) s;
      (void) tid;
      assert(false);
      return "";
    }

  /// \brief Returns an iterator over the successors of the parameter state.
  SuccIterator* succ(const State s, unsigned tid)
  {
    (void) s;
    (void) tid;
    assert(false);
    return nullptr;
  }

  /// \brief Allocation and deallocation of iterator is costly. This
  /// method allows to reuse old iterators.
  void recycle(SuccIterator* it, unsigned tid)
  {
    (void) it;
    (void) tid;
    assert(false);
  }

  /// \brief This method allow to deallocate a given state.
  const std::vector<std::string> ap()
  {
    assert(false);
    return {};
  }


  class PT_iterator
  {
  public:
    void next()
    {
      assert(false);
    }

    void done()
    {
      assert(false);
    }

    State state()
    {
      assert(false);
      return 0; // FIXME
    }

    spot::cube condition()
      {
	assert(false);
	return nullptr;
      }
  };

  // ------------------------------------------------------------------------
  // End matching the  Spot Kripkecube "interface"
  // ------------------------------------------------------------------------


};
