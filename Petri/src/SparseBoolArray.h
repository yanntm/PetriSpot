#ifndef SPARSEBOOLARRAY_H_
#define SPARSEBOOLARRAY_H_

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
 */

#include <vector>
#include <algorithm>
#include <iostream>

/**
 * SparseBoolArrays map integers to integers.  Unlike a normal array of integers,
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
class SparseBoolArray {
		int * mKeys;
    	int mSize;
		int mCap;
public:
	/**
     * Creates a new SparseBoolArray containing no mappings.
     */
    SparseBoolArray() : SparseBoolArray(10) {
    }
    /**
     * Creates a new SparseBoolArray containing no mappings that will not
     * require any additional memory allocation to store the specified
     * number of mappings.  If you supply an initial capacity of 0, the
     * sparse array will be initialized with a light-weight representation
     * not requiring any additional array allocations.
     */
    SparseBoolArray(int initialCapacity) {
		mCap = initialCapacity;
		mKeys = new int[mCap];
		mSize = 0;
    }
    /** 
     * Convert a classic std::vector<bool> to a sparse representation.
     * @param marks
     */
    SparseBoolArray(const std::vector<bool> & marks) {
	// compute and set correct capacity
		mCap = count_if(marks.begin(), marks.end(), [](const int & e) { return e != 0; });
		mKeys = new int[mCap];
		mSize = 0;
	   	for (int  i = 0, e = marks.size() ; i < e ; i++) {
	    	bool v = marks.at(i);
	    	if (v) {
		    	append(i, v);    			
	    	}
	   	}    	
    }

	~SparseBoolArray() {
		delete mKeys;
	}

    std::vector<bool> toList (int size) const {
    	std::vector<bool> res(size);
    	int  j = 0;
    	for (int i=0; i < size ; i++ ) {
    		if (j < this->size() && keyAt(j)==i) {
    			res.push_back(true);
    			++j;
    		} else {
    			res.push_back(false);
    		}    		
    	}
    	return res;
    }

	SparseBoolArray & operator=(const SparseBoolArray & source) {
		if (this != &source) {
			delete[] mKeys;
			mSize = mCap = 0;
			initFrom(*this, source);
		}
		return *this;
	}

private :
	static void initFrom(SparseBoolArray & a, const SparseBoolArray & o) {
		a.mCap = o.mCap;
		a.mSize = o.mSize;
		a.mKeys = new int[a.mCap];

		memcpy(a.mKeys, o.mKeys, a.mSize * sizeof(int));
	}

public :
	SparseBoolArray(const SparseBoolArray & o) {
		initFrom(*this,o);
	}
    
    
    /**
     * Gets the int mapped from the specified key, or <code>0</code>
     * if no such mapping has been made.
     */
    bool get(int key) const {
       	return get(key, false);
    }
    /**
     * Gets the int mapped from the specified key, or the specified value
     * if no such mapping has been made.
     */
    bool get(int key, bool valueIfKeyNotFound) const {
    	int i = binarySearch(mKeys, mSize, key);
       	if (i < 0) {
           		return valueIfKeyNotFound;
       	} else {
           		return true;
       	}
    }
    /**
     * Removes the mapping from the specified key, if there was any.
     */
    void remove(int key) {
       	int i = binarySearch(mKeys, mSize, key);
       	if (i >= 0) {
           		removeAt(i);
       	}
    }
    /**
     * Removes the mapping at the given index.
     */
    void removeAt(int index) {
       	std::copy(mKeys+index+1, mKeys+mSize, mKeys+index);  
       	mSize--;
    }
    /**
     * Adds a mapping from the specified key to the specified value,
     * replacing the previous mapping from the specified key if there
     * was one.
     */
    void put(int key, bool v) {    	
       	int i = binarySearch(mKeys, mSize, key);
       	if (i >= 0) {
       		if (v) {
       			return;
       		} else {
       			removeAt(i);
       		}
       	} else if (v) {
           		i = ~i;
           		auto p = insert(mKeys, mCap, mSize, i, key);
				mKeys = p.first;    
				mCap = p.second;   
           		mSize++;
       	}
    }
    /**
     * Returns the number of key-value mappings that this SparseBoolArray
     * currently stores.
     */
    int size() const {
       	return mSize;
    }
	/**
   	 * Given an index in the range <code>0...size()-1</code>, returns
   	 * the key from the <code>index</code>th key-value mapping that this
   	 * SparseBoolArray stores.
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
   	 * Returns the index for which {@link #keyAt} would return the
   	 * specified key, or a negative number if the specified
   	 * key is not mapped.
   	 */
   	int indexOfKey(int key) {
       	return binarySearch(mKeys, mSize, key);
   	}
   	/**
   	 * Removes all key-value mappings from this SparseBoolArray.
   	 */
   	void clear() {
       	mSize = 0;
   	}
   	/**
   	 * Puts a key/value pair into the array, optimizing for the case where
   	 * the key is greater than all existing keys in the array.
    	*/
   	void append(int key, bool v) {
   		if (! v)
   			return;
       	if (mSize != 0 && key <= mKeys[mSize - 1]) {
          	put(key, v);
        	return;
       	}
       	auto p = append(mKeys, mCap, mSize, key);
		mKeys = p.first;
		mCap = p.second;
      	mSize++;
    }
    /**
     * Provides a copy of keys.
     *
     * @hide
     * */
    int * copyKeys() {
       	if (size() == 0) {
           	return nullptr;
       	}
       	int* copiedKeys = new int[mSize];
        std::copy(mKeys, mKeys + mSize, copiedKeys);
        return copiedKeys;
    }
    /**
     * Provides direct access to keys, client should not modify.
     * @return an array of sorted integers corresponding to true entries of this BoolArray
     */
	int * refKeys() {
    	if (size() == 0) {
     		return new int[0];
    	}
    	return mKeys;
	}
    
    size_t hash() const {
		size_t prime = 1409;
		size_t result = 1;
		result = prime * result + hashCode(mKeys,mSize);
		result = prime * result + mSize;
		return result;
	}

    bool operator==(const SparseBoolArray & other) const {
		if (mSize != other.mSize)
			return false;

		if (!equalsRange(mKeys,other.mKeys,mSize))
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
     		return ;
       	}
     	os << '{';
    	for (int i=0; i<mSize; i++) {
    		if (i > 0) {
         		os <<  (", ");
        	}
       		int key = keyAt(i);
      		os << key;
    	}
       	os << '}';
       	return ;
    }

    friend std::ostream& operator<< (std::ostream & os, const SparseBoolArray & a) {
		a.print(os);
		return os;
	}
    
	void clear(int j) {
		put (j,false);
	}
	void set(int j) {
		put (j,true);
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

	static SparseBoolArray unionOperation(const SparseBoolArray& a, const SparseBoolArray& b) {
		int inter = 0;
		// step 1 evaluate size	
		int i=0;
		int j=0;
		for (int ie=a.size(), je=b.size() ; i < ie && j < je ; ) {
			if (a.mKeys[i] > b.mKeys[j]) {
				j++;
			} else if (a.mKeys[i] < b.mKeys[j]) {
				i++;
			} else {
				i++;
				j++;
				inter++;
			}
		}
		int resSize = a.size() + b.size() - inter;
		SparseBoolArray res(resSize);
		// step 2 assign
		int cur=0;
		i=0;
		j=0;
		for (int ie=a.size(), je=b.size() ; i < ie && j < je ; ) {
			if (a.mKeys[i] > b.mKeys[j]) {
				res.mKeys[cur++]=b.mKeys[j];
				j++;
			} else if (a.mKeys[i] < b.mKeys[j]) {
				res.mKeys[cur++]=a.mKeys[i];
				i++;
			} else {
				res.mKeys[cur++]=b.mKeys[j];
				i++;
				j++;
			}
		}
		// add remaining elements if any
		for (int ie=a.size(); i < ie ; i++) {
			res.mKeys[cur++]=a.mKeys[i];
		}
		for (int je=b.size(); j < je ; j++) {
			res.mKeys[cur++]=b.mKeys[j];
		}
		res.mSize = cur;
		return res;
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
		return currentSize <= 5 ? 10 : currentSize * 2;
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

	static int hashCode(const int * const a, int sz) {
		if (a == nullptr)
			return 0;

		int result = 1;
		for (int i=0; i < sz ; i++) {
			result = 31 * result + a[i];
		}

		return result;
	}

};


#endif /* SPARSEBOOLARRAY_H_ */
