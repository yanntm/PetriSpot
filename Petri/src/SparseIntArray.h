#ifndef SPARSEINTARRAY_H_
#define SPARSEINTARRAY_H_

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

#include <vector>
#include <utility>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <iostream>
#include <cmath>
#include <limits>

/**
 * SparseIntArrays map integers to integers.  Unlike a normal array of integers,
 * there can be gaps in the indices.  It is intended to be more memory efficient
 * than using a HashMap to map Integers to Integers, both because it avoids
 * auto-boxing keys and values and its data structure doesn't rely on an extra entry object
 * for each mapping.
 *
 * <p>Note that this container keeps its mappings in an array data structure,
 * using a binary search to find keys.  The implementation is not intended to be appropriate for
 * data structures
 * that may contain large numbers of items.  It is generally slower than a traditional
 * HashMap, since lookups require a binary search and adds and removes require inserting
 * and deleting entries in the array.  For containers holding up to hundreds of items,
 * the performance difference is not significant, less than 50%.</p>
 *
 * <p>It is possible to iterate over the items in this container using
 * {@link #keyAt(int)} and {@link #valueAt(int)}. Iterating over the keys using
 * <code>keyAt(int)</code> with ascending values of the index will return the
 * keys in ascending order, or the values corresponding to the keys in ascending
 * order in the case of <code>valueAt(int)</code>.</p>
 */
class SparseIntArray {
	int * mKeys;
	int * mValues;
	int mSize;
	int mCap;
public :
	/**
	 * Default constructor
	*/
	SparseIntArray() : mCap(8), mSize(0) {
        mKeys = new int[mCap];
        mValues = new int[mCap];
    }
	/**
	 * Creates a new SparseIntArray containing no mappings that will not
	 * require any additional memory allocation to store the specified
	 * number of mappings.  If you supply an initial capacity of 0, the
	 * sparse array will be initialized with a light-weight representation
	 * not requiring any additional array allocations.
	 */
	SparseIntArray(int initialCapacity=8) {
		mCap = initialCapacity;
		mKeys = new int[mCap];
		mValues = new int[mCap];
		mSize = 0;
	}
	/**
	 * Convert a classic vector<int> to a sparse representation.
	 * @param marks
	 */
	SparseIntArray(const std::vector<int> & marks) {
		// compute and set correct capacity
		mCap = count_if (marks.begin(), marks.end(), [](const int &e) { return e != 0; });
		mKeys = new int[mCap];
		mValues = new int[mCap];
		mSize = 0;
		for (int  i = 0, e = marks.size() ; i < e ; i++) {
			int v = marks.at(i);
			if (v != 0) {
				append(i, v);
			}
		}
	}

	~SparseIntArray () {
		delete [] mKeys;
		delete [] mValues;
	}

	std::vector<int> toList (int size) const {
		std::vector<int> res ;
		res.reserve(size);
		int  j = 0;
		for (int i=0; i < size ; i++ ) {
			if (j < this->size() && keyAt(j)==i) {
				res.push_back(valueAt(j));
				++j;
			} else {
				res.push_back(0);
			}
		}
		return res;
	}

//	/**
//	 * Acts like an std::move in c++ : steal the content of the argument.
//	 * This method updates the content of "this".
//	 * @param source an object that is invalidated by this operation.
//	 */
//	SparseIntArray (SparseIntArray && source) {
//		mKeys = source.mKeys;
//		source.mKeys = nullptr;
//		mValues = source.mValues;
//		source.mValues = nullptr;
//		mSize = source.mSize;
//		mCap = source.mCap;
//	}


	SparseIntArray & operator=(const SparseIntArray & source) {
		if (this != &source) {
			delete[] mKeys;
			delete[] mValues;
			mSize = mCap = 0;
			initFrom(*this, source);
		}
		return *this;
	}


private :
	static void initFrom(SparseIntArray & a, const SparseIntArray & o) {
		a.mCap = o.mCap;
		a.mSize = o.mSize;
		a.mKeys = new int[a.mCap];
		a.mValues = new int[a.mCap];;

		memcpy(a.mKeys, o.mKeys, a.mSize * sizeof(int));
		memcpy(a.mValues, o.mValues, a.mSize * sizeof(int));
	}

public :
	SparseIntArray (const SparseIntArray & o) {
		initFrom(*this,o);
	}
	/**
	 * Gets the int mapped from the specified key, or <code>0</code>
	 * if no such mapping has been made.
	 */
	int get(int key) const {
		return get(key, 0);
	}
	/**
	 * Gets the int mapped from the specified key, or the specified value
	 * if no such mapping has been made.
	 */
	int get(int key, int valueIfKeyNotFound) const {
		int i = binarySearch(mKeys, mSize, key);
		if (i < 0) {
			return valueIfKeyNotFound;
		} else {
			return mValues[i];
		}
	}
	/**
	 * Removes the mapping from the specified key, if there was any.
	 */
	int del(int key) {
		int i = binarySearch(mKeys, mSize, key);
		if (i >= 0) {
			removeAt(i);
		}
		return i;
	}
	/**
	 * Removes the mapping at the given index.
	 */
	void removeAt(int index) {
		std::copy(mKeys+index+1, mKeys+mSize, mKeys+index);
		std::copy(mValues+index+1, mValues+mSize, mValues+index);
		mSize--;
	}
	/**
	 * Adds a mapping from the specified key to the specified value,
	 * replacing the previous mapping from the specified key if there
	 * was one.
	 */
	void put(int key, int value) {
		int i = binarySearch(mKeys, mSize, key);
		if (value==0) {
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
				mValues = insert(mValues, mCap,  mSize, i, value).first;
				mCap = p.second;
				mSize++;
			}
		}
	}
	/**
	 * Returns the number of key-value mappings that this SparseIntArray
	 * currently stores.
	 */
	int size() const {
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
	int keyAt(int index) const {
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
	int valueAt(int index) const {
		return mValues[index];
	}
	/**
	 * Directly set the value at a particular index.
	 * @hide
	 */
	void setValueAt(int index, int value) {
		mValues[index] = value;
	}
	/**
	 * Returns the index for which {@link #keyAt} would return the
	 * specified key, or a negative number if the specified
	 * key is not mapped.
	 */
	int indexOfKey(int key) const {
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
	int indexOfValue(int value) const {
		for (int i = 0; i < mSize; i++)
			if (mValues[i] == value)
				return i;
		return -1;
	}
	/**
	 * Removes all key-value mappings from this SparseIntArray.
	 */
	void clear() {
		mSize = 0;
	}
	/**
	 * Puts a key/value pair into the array, optimizing for the case where
	 * the key is greater than all existing keys in the array.
	 */
	void append(int key, int value) {
		if (value == 0)
			return;
		if (mSize != 0 && key <= mKeys[mSize - 1]) {
			put(key, value);
			return;
		}
		auto p = append(mKeys, mCap, mSize, key);
		mKeys = p.first;
		mValues = append(mValues, mCap, mSize, value).first;
		mCap = p.second;
		mSize++;
	}

	size_t hash() const {
		size_t prime = 31;
		size_t result = 1;
		result = prime * result + hashCode(mKeys,mValues,mSize);
		result = prime * result + mSize;
		return result;
	}
	bool operator==(const SparseIntArray & other) const {
		if (mSize != other.mSize)
			return false;

		if (!equalsRange(mKeys,other.mKeys,mSize))
			return false;
		if (!equalsRange(mValues, other.mValues, mSize))
			return false;
		return true;
	}

	/**
	 * {@inheritDoc}
	 *
	 * <p>This implementation composes a string by iterating over its mappings.
	 */

	void print(std::ostream & os) const {
		if (size() <= 0) {
			os << "{}";
			return;
		}
		os <<  '{';
		for (int i=0; i<mSize; i++) {
			if (i > 0) {
				os <<  (", ");
			}
			int key = keyAt(i);
			os << key << "=" << valueAt(i);
		}
		os << '}';
		return ;
	}

	/**
	 * Returns a vector with : alpha * ta + beta * tb in each cell.
	 * @return
	 */
	static SparseIntArray sumProd(int alpha, const SparseIntArray & ta, int beta, const SparseIntArray & tb) {
		return sumProd(alpha, ta, beta, tb, -1);
	}
	static SparseIntArray sumProd(int alpha, const SparseIntArray & ta, int beta, const SparseIntArray & tb, int except) {
		SparseIntArray flow (std::max(ta.size(), tb.size()));

		int i = 0;
		int j = 0;
		while (i < ta.size() || j < tb.size()) {
			int ki = i==ta.size() ? std::numeric_limits<int>::max() : ta.keyAt(i);
			int kj = j==tb.size() ? std::numeric_limits<int>::max() : tb.keyAt(j);
			if (ki == kj) {
				int val = alpha * ta.valueAt(i)+ beta* tb.valueAt(j);
				if (val != 0 && ki != except) {
					flow.append(ki, val);
				}
				i++;
				j++;
			} else if (ki < kj) {
				int val = alpha * ta.valueAt(i);
				if (val != 0 && ki != except) flow.append(ki, val);
				i++;
			} else if (kj < ki) {
				int val = beta * tb.valueAt(j);
				if (val != 0 && kj != except) flow.append(kj, val);
				j++;
			}
		}

		return flow;
	}

	static bool keysIntersect(const SparseIntArray & s1, const SparseIntArray & s2) {
		if (s1.size() == 0 || s2.size() == 0) {
			return true;
		}

		for (int j = 0, i = 0 , ss1 =  s1.size() , ss2 = s2.size() ; i < ss1 && j < ss2 ; ) {
			int sk1 = s1.keyAt(i);
			int sk2 = s2.keyAt(j);
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
	 * Returns true iff. for all i, s1[i] >= s2[i].
	 * @param s1 the greater array
	 * @param s2 the smaller/lower valued array
	 * @return true iff for all i, s1[i] >= s2[i].
	 */
	static bool greaterOrEqual(const SparseIntArray & s1, const SparseIntArray & s2) {
		if (s1.size() < s2.size()) {
			return false;
		}
		if (s2.size() == 0) {
			return true;
		}


		for (int j = 0, i = 0 , ss1 =  s1.size() , ss2 = s2.size() ; i < ss1 && j < ss2 ; ) {
			int sk1 = s1.keyAt(i);
			int sk2 = s2.keyAt(j);
			if (sk1 == sk2) {
				if (s1.valueAt(i) < s2.valueAt(j)) {
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
				int ii = binarySearch(s1.mKeys, sk2, i+1, s1.mSize-1);
				if (ii < 0) {
					return false;
				}
				i = ii;
				if (ss1 - i  < ss2 - j) {
					return false;
				}
			}
		}
		return true;
	}

	static int manhattanDistance (const SparseIntArray& ta, const SparseIntArray& tb) {
		int dist = 0;

		int i = 0;
		int j = 0;
		while (i < ta.size() || j < tb.size()) {
			int ki = i==ta.size() ? std::numeric_limits<int>::max() : ta.keyAt(i);
			int kj = j==tb.size() ? std::numeric_limits<int>::max() : tb.keyAt(j);
			if (ki == kj) {
				dist += abs(ta.valueAt(i) - tb.valueAt(j));
				i++;
				j++;
			} else if (ki < kj) {
				int val = ta.valueAt(i);
				dist += abs(val);
				i++;
			} else if (kj < ki) {
				int val = tb.valueAt(j);
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
	void deleteAndShift(int i) {
		if (mSize==0 || i > mKeys[mSize-1]) {
			return;
		}
		int k;
		for (k= mSize-1 ; k>=0 && mKeys[k]>i ; k--) {
			mKeys[k]--;
		}
		if (k >= 0 && mKeys[k]==i) {
			removeAt(k);
		}
	}

	/**
	 * More efficient one pass version for removing many at once.
	 * @param todel a list of decreasing sorted indexes.
	 */
	void deleteAndShift(const std::vector<int> & todel) {
		if (mSize==0 || todel.empty() || todel.back() > mKeys[mSize-1]) {
			return;
		}
		// start from rightmost key
		int k = mSize-1 ;
		for (int cur=0, e=todel.size() ; cur < e ; cur++) {
			int val = todel.at(cur);
			for ( ; k>=0 && mKeys[k]> val ; k--) {
				mKeys[k]-= todel.size() - cur;
			}
			if (k >= 0 && mKeys[k]==val) {
				removeAt(k);
				--k;
			}
		}
	}

private :
	static int binarySearch(const int * const array, int sz, int value) {
		int lo = 0;
		int hi = sz - 1;

		return binarySearch(array, value, lo, hi);
	}
	// This is Arrays.binarySearch(), but doesn't do any argument validation.
	static int binarySearch(const int * const array, int value, int lo, int hi) {
		while (lo <= hi) {
			int mid = (lo + hi) >> 1;
			int midVal = array[mid];
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
	static int growSize(int currentSize) {
		return currentSize <= 4 ? 8 : currentSize * 2;
	}
	static std::pair<int*,int> append(int* array, int aCap, int currentSize, int element) {
		assert (currentSize <= aCap) ;
		int nCap = aCap;
		if (currentSize + 1 > aCap) {
			nCap = growSize(currentSize);
			int * newArray = new int [nCap];
			memcpy(newArray, array, currentSize * sizeof(int));
			delete [] array;
			array = newArray;
		}
		array[currentSize] = element;
		return {array,nCap};
	}

	/**
	 * Primitive int version of {@link #insert(Object[], int, int, Object)}.
	 */
	static std::pair<int*,int> insert(int* array, int aCap, int currentSize, int index, int element) {
		assert (currentSize <= aCap);
		if (currentSize + 1 <= aCap) {
			memmove(array+index+1, array+index, (currentSize - index)*sizeof(int));
			array[index] = element;
			return {array,aCap};
		}
		int nCap = growSize(currentSize);
		int * newArray = new int [nCap];
		memcpy(newArray, array, index * sizeof(int));
		newArray[index] = element;
		memcpy(newArray+index+1, array+index, (aCap -index) * sizeof(int));
		delete [] array;
		return {newArray,nCap};
	}

	static bool equalsRange(const int* const a, const int* const b, int s) {
		return std::equal(a, a +s, b);
	}

	static int hashCode(const int * const a, const int * const b, int sz) {
		if (a == nullptr || b==nullptr)
			return 0;

		int result = 1;
		for (int i=0; i < sz ; i++) {
			result = 31 * result + a[i];
			result = 31 * result + b[i];
		}

		return result;
	}

};



#endif /* SPARSEINTARRAY_H_ */
