#ifndef SPARSEARRAY_T_H_
#define SPARSEARRAY_T_H_

/*
 * Copyright (C) 2006 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * package android.util;
 */
#include <functional>
#include <vector>
#include <utility>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <iostream>
#include <cmath>
#include <limits>
#include <concepts>

/**
 * @brief SparseArray maps integers to values of type T. Unlike a normal array,
 * there can be gaps in the indices. It is intended to be more memory efficient
 * than using a HashMap to map integers to values, as it avoids auto-boxing keys
 * and values and does not rely on an extra entry object for each mapping.
 *
 * Note that this container keeps its mappings in an array data structure,
 * using a binary search to find keys. The implementation is not intended for
 * data structures that may contain large numbers of items. It is generally
 * slower than a traditional HashMap, since lookups require a binary search and
 * additions and removals require inserting and deleting entries in the array.
 * For containers holding up to hundreds of items, the performance difference
 * is not significant, typically less than 50%.
 *
 * It is possible to iterate over the items in this container using
 * keyAt(size_t) and valueAt(size_t). Iterating over the keys using
 * keyAt(size_t) with ascending values of the index will return the keys in
 * ascending order, and the values corresponding to the keys in ascending
 * order in the case of valueAt(size_t).
 *
 * @tparam T The type of the values to be stored in the array.
 */
template<std::integral T>
class SparseArray {
    size_t *mKeys;
    T *mValues;
    size_t mSize;
    size_t mCap;

public:
    /**
     * Creates a new SparseIntArray containing no mappings that will not
     * require any additional memory allocation to store the specified
     * number of mappings. If you supply an initial capacity of 0, the
     * sparse array will be initialized with a light-weight representation
     * not requiring any additional array allocations.
     */
    SparseArray(size_t initialCapacity = 0) {
        mCap = initialCapacity;
        mSize = 0;
        if (mCap == 0) {
            mKeys = nullptr;
            mValues = nullptr;
        } else {
            mKeys = new size_t[mCap];
            mValues = new T[mCap];
        }
    }

    /**
     * Convert a classic vector<U> to a sparse representation.
     * NB : U should be comparable by != to 0 and copy assignable to T.
     * @param marks Vector of type U elements.
     */
    template<typename U>
    SparseArray(const std::vector<U> &marks) {
        mCap = std::count_if(marks.begin(), marks.end(),
                             [](const U &e) { return e != 0; });
        mSize = 0;
        if (mCap == 0) {
            mKeys = nullptr;
            mValues = nullptr;
        } else {
            mKeys = new size_t[mCap];
            mValues = new T[mCap];
            for (size_t i = 0, e = marks.size(); i < e; ++i) {
                U v = marks.at(i);
                if (v != 0) {
                    append(i, static_cast<T>(v));
                }
            }
        }
    }

    SparseArray(SparseArray&& other) noexcept
        : mKeys(nullptr), mValues(nullptr), mSize(0), mCap(0) {
        mKeys = other.mKeys;
        mValues = other.mValues;
        mSize = other.mSize;
        mCap = other.mCap;

        other.mKeys = nullptr;
        other.mValues = nullptr;
        other.mSize = 0;
        other.mCap = 0;
    }

    SparseArray& operator=(SparseArray&& other) noexcept {
        if (this != &other) {
            delete[] mKeys;
            delete[] mValues;

            mKeys = other.mKeys;
            mValues = other.mValues;
            mSize = other.mSize;
            mCap = other.mCap;

            other.mKeys = nullptr;
            other.mValues = nullptr;
            other.mSize = 0;
            other.mCap = 0;
        }
        return *this;
    }

    ~SparseArray() {
        delete[] mKeys;
        delete[] mValues;
    }

    std::vector<T> toList(size_t size) const {
        std::vector<T> res;
        res.reserve(size);
        size_t j = 0;
        for (size_t i = 0; i < size; i++) {
            if (j < this->size() && keyAt(j) == i) {
                res.push_back(valueAt(j));
                ++j;
            } else {
                res.push_back(0);
            }
        }
        return res;
    }

    void scalarProduct(T factor) {
        for (size_t i = 0; i < mSize; i++) {
            mValues[i] *= factor;
        }
    }

    void scalarDiv(T factor) {
        for (size_t i = 0; i < mSize; i++) {
            mValues[i] /= factor;
        }
    }

    /**
     * Acts like an std::move in C++: steal the content of the argument.
     * This method updates the content of "this".
     * @param source an object that is invalidated by this operation.
     */
    void move(SparseArray &source) {
        delete[] mKeys;
        mKeys = source.mKeys;
        source.mKeys = nullptr;
        delete[] mValues;
        mValues = source.mValues;
        source.mValues = nullptr;
        mSize = source.mSize;
        mCap = source.mCap;
    }

    SparseArray& operator=(const SparseArray &source) {
        if (this == &source) {
            return *this;
        }

        // Rule 1: If source is empty, clear to lightweight state
        if (source.mSize == 0) {
            delete[] mKeys;
            delete[] mValues;
            mKeys = nullptr;
            mValues = nullptr;
            mSize = 0;
            mCap = 0;
            return *this;
        }

        // Rule 2: If mCap is sufficient reuse it
        // Commented out : if excessively large, realloc
        if (mCap >= source.mSize /* && mCap < 2 * source.mSize*/) {
            memcpy(mKeys, source.mKeys, source.mSize * sizeof(size_t));
            memcpy(mValues, source.mValues, source.mSize * sizeof(T));
            mSize = source.mSize;
            return *this;
        }

        // Rule 3: Deallocate, allocate to source.mSize, copy data
        delete[] mKeys;
        delete[] mValues;
        mCap = source.mSize;
        mKeys = new size_t[mCap];
        mValues = new T[mCap];
        memcpy(mKeys, source.mKeys, source.mSize * sizeof(size_t));
        memcpy(mValues, source.mValues, source.mSize * sizeof(T));
        mSize = source.mSize;
        return *this;
    }

private:
    static void initFrom(SparseArray &a, const SparseArray &o) {
        a.mCap = o.mCap;
        a.mSize = o.mSize;
        if (a.mCap == 0) {
            a.mKeys = nullptr;
            a.mValues = nullptr;
        } else {
            a.mKeys = new size_t[a.mCap];
            a.mValues = new T[a.mCap];
            memcpy(a.mKeys, o.mKeys, a.mSize * sizeof(size_t));
            memcpy(a.mValues, o.mValues, a.mSize * sizeof(T));
        }
    }

public:
    SparseArray(const SparseArray &o) {
        initFrom(*this, o);
    }

    /**
     * Gets the value mapped from the specified key, or <code>0</code>
     * if no such mapping has been made.
     */
    T get(size_t key) const {
        return get(key, 0);
    }

    /**
     * Gets the value mapped from the specified key, or the specified value
     * if no such mapping has been made.
     */
    T get(size_t key, const T &valueIfKeyNotFound) const {
        ssize_t i = binarySearch(mKeys, mSize, key);
        if (i < 0) {
            return valueIfKeyNotFound;
        } else {
            return mValues[i];
        }
    }

    /**
     * Removes the mapping from the specified key, if there was any.
     */
    ssize_t del(size_t key) {
        ssize_t i = binarySearch(mKeys, mSize, key);
        if (i >= 0) {
            removeAt(i);
        }
        return i;
    }

    /**
     * Removes the mapping at the given index.
     */
    void removeAt(size_t index) {
        std::copy(mKeys + index + 1, mKeys + mSize, mKeys + index);
        std::copy(mValues + index + 1, mValues + mSize, mValues + index);
        mSize--;
    }

    /**
     * Adds a mapping from the specified key to the specified value,
     * replacing the previous mapping from the specified key if there
     * was one.
     */
    void put(size_t key, T value) {
        ssize_t i = binarySearch(mKeys, mSize, key);
        if (value == 0) {
            if (i >= 0) {
                removeAt(i);
            }
        } else {
            if (i >= 0) {
                mValues[i] = value;
            } else {
                i = ~i;
                auto p = insert(mKeys, mCap, mSize, i, key);
                mKeys = p.first;
                mValues = insert(mValues, mCap, mSize, i, value).first;
                mCap = p.second;
                mSize++;
            }
        }
    }

    /**
     * Returns the number of key-value mappings that this SparseIntArray
     * currently stores.
     */
    size_t size() const {
        return mSize;
    }

    /**
     * Given an index in the range <code>0...size()-1</code>, returns
     * the key from the <code>index</code>th key-value mapping that this
     * SparseIntArray stores.
     *
     * <p>The keys corresponding to indices in ascending order are guaranteed to
     * be in ascending order, e.g., <code>keyAt(0)</code> will return the
     * smallest key and <code>keyAt(size()-1)</code> will return the largest
     * key.</p>
     */
    size_t keyAt(size_t index) const {
        assert(index < mSize);
        return mKeys[index];
    }

    /**
     * Given an index in the range <code>0...size()-1</code>, returns
     * the value from the <code>index</code>th key-value mapping that this
     * SparseIntArray stores.
     *
     * <p>The values corresponding to indices in ascending order are guaranteed
     * to be associated with keys in ascending order, e.g.,
     * <code>valueAt(0)</code> will return the value associated with the
     * smallest key and <code>valueAt(size()-1)</code> will return the value
     * associated with the largest key.</p>
     */
    T valueAt(size_t index) const {
        assert(index < mSize);
        return mValues[index];
    }

    /**
     * Directly set the value at a particular index.
     * @hide
     */
    void setValueAt(size_t index, T value) {
        assert(index < mSize);
        mValues[index] = value;
        // If setting to 0 and this was the last non-zero element, clear
        if (value == 0 && mSize == 1) {
            clear();
        }
    }

    /**
     * Returns the index for which {@link #keyAt} would return the
     * specified key, or a negative number if the specified
     * key is not mapped.
     */
    ssize_t indexOfKey(size_t key) const {
        return binarySearch(mKeys, mSize, key);
    }

    /**
     * Returns an index for which {@link #valueAt} would return the
     * specified key, or a negative number if no keys map to the
     * specified value.
     * Beware that this is a linear search, unlike lookups by key,
     * and that multiple keys can map to the same value and this will
     * find only one of them.
     */
    int indexOfValue(T value) const {
        for (size_t i = 0; i < mSize; i++) {
            if (mValues[i] == value) return i;
        }
        return -1;
    }

    /**
     * Removes all key-value mappings from this SparseIntArray.
     */
    void clear() {
        delete[] mKeys;
        delete[] mValues;
        mKeys = nullptr;
        mValues = nullptr;
        mSize = 0;
        mCap = 0;
    }

    // Reset size to 0, keep capacity for reuse
    void reset() {
        mSize = 0;
    }

    /**
     * Puts a key/value pair into the array, optimizing for the case where
     * the key is greater than all existing keys in the array.
     */
    void append(size_t key, T value) {
        if (value == 0) return;
        if (mSize != 0 && key <= mKeys[mSize - 1]) {
            put(key, value);
            return;
        }
        auto p = append(mKeys, mCap, mSize, key);
        mKeys = p.first;
        auto p2 = append(mValues, mCap, mSize, value);
        mValues = p2.first;
        mCap = p.second;
        mSize++;
    }

    size_t hash() const {
        size_t prime = 31;
        size_t result = 1;
        result = prime * result + hashCode(mKeys, mValues, mSize);
        result = prime * result + mSize;
        return result;
    }

    bool operator==(const SparseArray &other) const {
        if (mSize != other.mSize) return false;
        if (memcmp(mKeys, other.mKeys, mSize * sizeof(size_t)) != 0) {
            return false;
        }
        if (memcmp(mValues, other.mValues, mSize * sizeof(T)) != 0) {
            return false;
        }
        return true;
    }

    T maxVal() const {
        if (size() == 0) return 0;
        T max = mValues[0];
        for (size_t i = 1; i < mSize; i++) {
            if (std::abs(mValues[i]) > std::abs(max)) {
                max = mValues[i];
            }
        }
        return max;
    }

    /**
     * {@inheritDoc}
     *
     * <p>This implementation composes a string by iterating over its mappings.
     */
    void print(std::ostream &os) const {
        if (size() <= 0) {
            os << "{}";
            return;
        }
        os << '{';
        for (size_t i = 0; i < mSize; i++) {
            if (i > 0) {
                os << ", ";
            }
            size_t key = keyAt(i);
            os << key << "=" << valueAt(i);
        }
        os << '}';
    }

    friend std::ostream& operator<<(std::ostream &os, const SparseArray &a) {
        a.print(os);
        return os;
    }

    /**
     * Returns a vector with: alpha * ta + beta * tb in each cell.
     * @return
     */
    static SparseArray sumProd(int alpha, const SparseArray &ta, int beta,
                               const SparseArray &tb) {
        return sumProd(alpha, ta, beta, tb, -1);
    }

    static SparseArray sumProd(int alpha, const SparseArray &ta, int beta,
                               const SparseArray &tb, int except) {
        SparseArray flow(std::max(ta.size(), tb.size()));

        size_t i = 0;
        size_t j = 0;
        while (i < ta.size() || j < tb.size()) {
            size_t ki = i == ta.size() ? std::numeric_limits<size_t>::max() : ta.keyAt(i);
            size_t kj = j == tb.size() ? std::numeric_limits<size_t>::max() : tb.keyAt(j);
            if (ki == kj) {
                T val = alpha * ta.valueAt(i) + beta * tb.valueAt(j);
                if (val != 0 && (except > 0 ? ki != (size_t)except : 1)) {
                    flow.append(ki, val);
                }
                i++;
                j++;
            } else if (ki < kj) {
                T val = alpha * ta.valueAt(i);
                if (val != 0 && (except > 0 ? ki != (size_t)except : 1)) {
                    flow.append(ki, val);
                }
                i++;
            } else if (kj < ki) {
                T val = beta * tb.valueAt(j);
                if (val != 0 && (except > 0 ? kj != (size_t)except : 1)) {
                    flow.append(kj, val);
                }
                j++;
            }
        }
        return flow;
    }

    bool isPurePositive() const {
        for (size_t i = 0; i < mSize; i++) {
            if (mValues[i] < 0) return false;
        }
        return true;
    }

    static bool keysIntersect(const SparseArray &s1, const SparseArray &s2) {
        if (s1.size() == 0 || s2.size() == 0) {
            return true;
        }
        for (size_t j = 0, i = 0, ss1 = s1.size(), ss2 = s2.size(); i < ss1 && j < ss2;) {
            size_t sk1 = s1.keyAt(i);
            size_t sk2 = s2.keyAt(j);
            if (sk1 == sk2) {
                return true;
            } else if (sk1 > sk2) {
                j++;
            } else {
                i++;
            }
        }
        return false;
    }

    /**
     * Test whether all entries in s1 are greater or equal than the corresponding entry in s2.
     * Returns true iff for all i, s1[i] >= s2[i].
     * @param s1 the greater array
     * @param s2 the smaller/lower valued array
     * @return true iff for all i, s1[i] >= s2[i].
     */
    static bool greaterOrEqual (const SparseArray &s1, const SparseArray &s2)
    {
      if (s1.size () < s2.size ()) {
        return false;
      }
      if (s2.size () == 0) {
        return true;
      }

      for (size_t j = 0, i = 0, ss1 = s1.size (), ss2 = s2.size ();
          i < ss1 && j < ss2;) {
        size_t sk1 = s1.keyAt (i);
        size_t sk2 = s2.keyAt (j);
        if (sk1 == sk2) {
          if (s1.valueAt (i) < s2.valueAt (j)) {
            return false;
          } else {
            i++;
            j++;
          }
        } else if (sk1 > sk2) {
          // missing entries !
          return false;
        } else {
          // sk1 < sk2 : we must progress in s1
          // use a binary search for that
          ssize_t ii = binarySearch (s1.mKeys, sk2, i + 1, s1.mSize - 1);
          if (ii < 0) {
            return false;
          }
          i = ii;
          if (ss1 - i < ss2 - j) {
            return false;
          }
        }
      }
      return true;
    }

    /**
     * Computes the maximum number of times s2 can be subtracted from s1, assuming all entries are positive.
     * Returns the minimum integer ratio s1[i]/s2[i] over common keys, or 0 if s2 is not contained in s1.
     * @param s1 the greater array
     * @param s2 the smaller/lower valued array
     * @return the integer k such that s1 >= k * s2, or 0 if not fully contained
     */
    static ssize_t countContainsPos (const SparseArray<T> &s1, const SparseArray<T> &s2)
    {
      if (s1.size () < s2.size ()) {
        return 0;
      }
      if (s2.size () == 0) {
        return 0;
      }

      ssize_t minCount = std::numeric_limits<ssize_t>::max();
      for (ssize_t j = 0, i = 0, ss1 = s1.size (), ss2 = s2.size ();
          i < ss1 && j < ss2;) {
        ssize_t sk1 = s1.keyAt (i);
        ssize_t sk2 = s2.keyAt (j);
        if (sk1 == sk2) {
          T count = s1.valueAt (i) / s2.valueAt (j);
          if (count == 0) {
            return 0;
          } else {
            if (count < minCount) {
              minCount = count;
            }
            i++;
            j++;
          }
        } else if (sk1 > sk2) {
          // missing entries !
          return 0;
        } else {
          // sk1 < sk2 : we must progress in s1
          // use a binary search for that
          ssize_t ii = binarySearch (s1.mKeys, sk2, i + 1, s1.mSize - 1);
          if (ii < 0) {
            return 0;
          }
          i = ii;
          if (ss1 - i < ss2 - j) {
            return 0;
          }
        }
      }
      return minCount == std::numeric_limits<ssize_t>::max () ? 0 : minCount;
    }

    static size_t manhattanDistance(const SparseArray &ta, const SparseArray &tb) {
        size_t dist = 0;
        size_t i = 0;
        size_t j = 0;
        while (i < ta.size() || j < tb.size()) {
            size_t ki = i == ta.size() ? std::numeric_limits<size_t>::max() : ta.keyAt(i);
            size_t kj = j == tb.size() ? std::numeric_limits<size_t>::max() : tb.keyAt(j);
            if (ki == kj) {
                dist += abs(ta.valueAt(i) - tb.valueAt(j));
                i++;
                j++;
            } else if (ki < kj) {
                T val = ta.valueAt(i);
                dist += abs(val);
                i++;
            } else if (kj < ki) {
                T val = tb.valueAt(j);
                dist += abs(val);
                j++;
            }
        }
        return dist;
    }

    /**
     * Delete an element at index and shift elements to the right by one.
     * @param i
     */
    void deleteAndShift(size_t i) {
        if (mSize == 0 || i > mKeys[mSize - 1]) {
            return;
        }
        ssize_t k = static_cast<ssize_t>(mSize) - 1;
        for (; k >= 0 && mKeys[k] > i; k--) {
            mKeys[k]--;
        }
        if (k >= 0 && mKeys[k] == i) {
            removeAt(static_cast<size_t>(k));
        }
    }

    /**
     * More efficient one-pass version for removing many at once.
     * @param todel a list of decreasing sorted indexes.
     */
    void deleteAndShift(const std::vector<size_t> &todel) {
        if (mSize == 0 || todel.empty() || todel.back() > mKeys[mSize - 1]) {
            return;
        }
        ssize_t k = static_cast<ssize_t>(mSize) - 1;
        for (size_t cur = 0, e = todel.size(); cur < e; cur++) {
            size_t val = todel[cur];
            for (; k >= 0 && mKeys[k] > val; k--) {
                mKeys[k] -= (todel.size() - cur);
            }
            if (k >= 0 && mKeys[k] == val) {
                removeAt(static_cast<size_t>(k));
                --k;
            }
        }
    }

private:
    static ssize_t binarySearch(const size_t *const array, size_t sz, size_t value) {
        ssize_t lo = 0;
        ssize_t hi = sz - 1;
        return binarySearch(array, value, lo, hi);
    }

    static ssize_t binarySearch(const size_t *const array, size_t value,
                                       ssize_t lo, ssize_t hi) {
        while (lo <= hi) {
            ssize_t mid = (lo + hi) >> 1;
            size_t midVal = array[mid];
            if (midVal < value) {
                lo = mid + 1;
            } else if (midVal > value) {
                hi = mid - 1;
            } else {
                return mid;  // value found
            }
        }
        return ~lo;  // value not present
    }

    /**
     * Given the current size of an array, returns an ideal size to which the array should grow.
     * This is typically double the given size, but should not be relied upon to do so in the
     * future.
     */
    static size_t growSize(size_t currentSize) {
        return currentSize <= 4 ? 8 : currentSize * 2;
    }

    template<typename U>
    static std::pair<U*, size_t> append(U *array, size_t aCap, size_t currentSize, U element) {
        assert(currentSize <= aCap);
        size_t nCap = aCap;
        if (currentSize + 1 > aCap) {
            nCap = growSize(currentSize);
            U *newArray = new U[nCap];
            if (array) {
                memcpy(newArray, array, currentSize * sizeof(U));
                delete[] array;
            }
            array = newArray;
        }
        array[currentSize] = element;
        return {array, nCap};
    }

    template<typename U>
    static std::pair<U*, size_t> insert(U *array, size_t aCap, size_t currentSize,
                                        size_t index, U element) {
        assert(currentSize <= aCap);
        if (currentSize + 1 <= aCap) {
            memmove(array + index + 1, array + index, (currentSize - index) * sizeof(U));
            array[index] = element;
            return {array, aCap};
        }
        size_t nCap = growSize(currentSize);
        U *newArray = new U[nCap];
        if (array) {
            memcpy(newArray, array, index * sizeof(U));
            memcpy(newArray + index + 1, array + index, (aCap - index) * sizeof(U));
            delete[] array;
        }
        newArray[index] = element;
        return {newArray, nCap};
    }

    static size_t hashCode(const size_t *const a, const T *const b, size_t sz) {
        if (a == nullptr || b == nullptr) return 0;
        size_t result = 1;
        for (size_t i = 0; i < sz; i++) {
            result = 31 * result + a[i];
            result = 31 * result + b[i];
        }
        return result;
    }
};

// Specialize std::hash for SparseArray<T>
namespace std {
template<typename T>
struct hash<SparseArray<T>> {
    size_t operator()(const SparseArray<T> &obj) const {
        return obj.hash();
    }
};

template<typename T>
struct hash<SparseArray<T>*> {
    size_t operator()(const SparseArray<T> *ptr) const {
        return ptr ? ptr->hash() : 0;
    }
};

template<typename T>
struct equal_to<SparseArray<T>*> {
    bool operator()(const SparseArray<T> *lhs, const SparseArray<T> *rhs) const {
        if (lhs == rhs) return true;
        if (!lhs || !rhs) return false;
        return *lhs == *rhs;
    }
};

template<typename T>
struct hash<const SparseArray<T>*> {
    size_t operator()(const SparseArray<T> *ptr) const {
        return ptr ? ptr->hash() : 0;
    }
};

template<typename T>
struct equal_to<const SparseArray<T>*> {
    bool operator()(const SparseArray<T> *lhs, const SparseArray<T> *rhs) const {
        if (lhs == rhs) return true;
        if (!lhs || !rhs) return false;
        return *lhs == *rhs;
    }
};
}

typedef SparseArray<int> SparseIntArray;

#endif /* SPARSEARRAY_T_H_ */
