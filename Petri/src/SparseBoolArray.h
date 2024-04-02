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
class SparseBoolArray {
private:
    std::vector<int> mKeys;
    int mSize;
public:
    /**
     * Creates a new SparseIntArray containing no mappings.
     */
    SparseBoolArray() : SparseBoolArray(10) {
    }
    /**
     * Creates a new SparseIntArray containing no mappings that will not
     * require any additional memory allocation to store the specified
     * number of mappings.  If you supply an initial capacity of 0, the
     * sparse array will be initialized with a light-weight representation
     * not requiring any additional array allocations.
     */
    SparseBoolArray(int initialCapacity) : mKeys(std::vector<int>(initialCapacity)) , mSize(0) {
    }
    /** 
     * Convert a classic List<Int> to a sparse representation.
     * @param marks
     */
    SparseBoolArray(std::vector<bool> marks) : SparseBoolArray(std::count_if(marks.begin(), marks.end(), [](int e) { return e != 0; })) {
    	// compute and set correct capacity
    	for (int  i = 0, e = marks.size() ; i < e ; i++) {
    		bool v = marks[i];
    		if (v) {
    			append(i, v);    			
    		}
    	}    	
	}

    std::vector<bool> toList (int size) const {
    	std::vector<bool> res(size);
    	int  j = 0;
    	for (int i=0; i < size ; i++ ) {
    		if (j < SparseBoolArray::size() && keyAt(j)==i) {
    			res.push_back(true);
    			++j;
    		} else {
    			res.push_back(false);
    		}    		
    	}
    	return res;
    }
    
    SparseBoolArray clone() {
        SparseBoolArray clone;
        try {
            clone.mSize = mSize;
            clone.mKeys = mKeys;
        } catch (const char* e) {
		    std::cout << e << std::endl;
		    return 1;
	    }
        return clone;
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
        auto it = std::find(mKeys.begin(), mKeys.end(), key);
        if (it != mKeys.end()) {
            return true;
        } else {
            return valueIfKeyNotFound;
        }
    }
    /**
     * Removes the mapping from the specified key, if there was any.
     */
    void remove(int key) {
        auto it = std::find(mKeys.begin(), mKeys.end(), key);
        if (it != mKeys.end()) {
            removeAt(std::distance(mKeys.begin(), it));
        }
    }

// poursuivre ici

    /**
     * Removes the mapping at the given index.
     */
    void removeAt(int index) {
        System.arraycopy(mKeys, index + 1, mKeys, index, mSize - (index + 1));        
        mSize--;
    }
    /**
     * Adds a mapping from the specified key to the specified value,
     * replacing the previous mapping from the specified key if there
     * was one.
     */
    void put(int key, bool v) {    	
        auto it = std::find(mKeys.begin(), mKeys.end(), key);
        if (it != mKeys.end()) {
        	if (v) {
        		return;
        	} else {
        		removeAt(std::distance(mKeys.begin(), it));
        	}
        } else if (v) {
            i = ~i;
            mKeys = GrowingArrayUtils.insert(mKeys, mSize, i, key);            
            mSize++;
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
     * Returns the index for which {@link #keyAt} would return the
     * specified key, or a negative number if the specified
     * key is not mapped.
     */
    int indexOfKey(int key) {
        return std::distance(mKeys.begin(), std::find(mKeys.begin(), mKeys.end(), key));
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
        mKeys.push_back(key);
        mSize++;
    }
    /**
     * Provides a copy of keys.
     *
     * @hide
     * */
    public int[] copyKeys() {
        if (size() == 0) {
            return null;
        }
        return Arrays.copyOf(mKeys, size());
    }
    /**
     * Provides direct access to keys, client should not modify.
     * As side effect, trims the representation array if it's overlarge.
     * @return an array of sorted integers corresponding to true entries of this BoolArray
     */
	public int[] refKeys() {
        if (size() == 0) {
            return new int[0];
        }
        if (mKeys.length > size()) {
        	mKeys = Arrays.copyOf(mKeys, size()); 
        }
        return mKeys;
	}

    
    @Override
	public int hashCode() {
		final int prime = 1409;
		int result = 1;
		result = prime * result + ContainerHelpers.hashCode(mKeys,mSize);
		result = prime * result + mSize;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof SparseBoolArray))
			return false;
		SparseBoolArray other = (SparseBoolArray) obj;		
		if (mSize != other.mSize)
			return false;

		if (!equalsRange(mKeys,other.mKeys,mSize))
			return false;
		return true;
	}
	private boolean equalsRange(int[] a, int[] b, int s) {
		if (a==b) {
			return true;
		}
		for (int i=0; i< s; i++) {
			if (a[i] != b[i])
				return false;
		}
		return true;
	}
	/**
     * {@inheritDoc}
     *
     * <p>This implementation composes a string by iterating over its mappings.
     */
    @Override
    public String toString() {
        if (size() <= 0) {
            return "{}";
        }
        StringBuilder buffer = new StringBuilder(mSize * 28);
        buffer.append('{');
        for (int i=0; i<mSize; i++) {
            if (i > 0) {
                buffer.append(", ");
            }
            int key = keyAt(i);
            buffer.append(key);
        }
        buffer.append('}');
        return buffer.toString();
    }
    
	public void clear(int j) {
		put (j,false);
	}
	public void set(int j) {
		put (j,true);
	}
	
	 /**
     * Delete an element at index and shift elements to the right by one.
     * @param i
     */
	public void deleteAndShift(int i) {
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
	public static SparseBoolArray or(SparseBoolArray a, SparseBoolArray b) {
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
		SparseBoolArray res = new SparseBoolArray(resSize);
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


};


#endif /* SPARSEBOOLARRAY_H_ */