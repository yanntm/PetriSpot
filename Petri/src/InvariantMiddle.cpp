#include "InvariantMiddle.h"


std::mutex InvariantMiddle::lock;
MatrixCol InvariantMiddle::last;
std::unordered_set<SparseIntArray> InvariantMiddle::lastInv;