#ifndef INVARIANTCALCULATOR_H_
#define INVARIANTCALCULATOR_H_

/**
 * A calculator for invariants and testing if a net is covered by invariants.
 * Provides two differient algorithms for calculating the invariants. The first
 * algorithm is descripted in http://de.scribd.com/doc/49919842/Pn-ESTII (slide
 * 88) and the other is also based an the farkas algorithm and is descripted in
 * http://pipe2.sourceforge.net/documents/PIPE-Report.pdf (page 19) which is
 * based on the paper of D'Anna and Trigila "Concurrent system analysis using
 * Petri nets – an optimised algorithm for finding net invariants", Mario D'Anna
 * and Sebastiano Trigila, Computer Communications vol 11, no. 4 august 1988.
 *
 * @author Dennis-Michael Borde, Manuel Gieseking , Adapted to ITS-tools by Yann
 *         Thierry-Mieg, 2017.
 *
 */

#include <string>
#include <iostream>
#include <vector>
#include "MatrixCol.h"

class InvariantCalculator {

	static const bool DEBUG = false;

	/**
	 * Hidden constructor
	 */
private:	
	InvariantCalculator() {
	}

	/**
	 * A class for holding the sets P+ = {j | c_hj &gt; 0} and P- = {j | c_hj less
	 * 0} for a given row.
	 */
	class PpPm {

	public:
		// The row
		const int row;
		// P+ set
		const SparseBoolArray pPlus;
		// P- set
		const SparseBoolArray pMinus;

		/**
		 * initially empty.
		 *
		 * @param row
		 */
		PpPm(int row) : row(row), pPlus(new SparseBoolArray()), pMinus(new SparseBoolArray()) {
		}

		void setValue(int j, int val) {
			if (val == 0) {
				pMinus.clear(j);
				pPlus.clear(j);
			} else if (val < 0) {
				pMinus.set(j);
				pPlus.clear(j);
			} else {
				pMinus.clear(j);
				pPlus.set(j);
			}
		}

		friend std::ostream& operator<<(std::ostream& os, const PpPm& obj) {
			os << "PpPm [row=" << obj.row << ", pPlus=" << obj.pPlus << ", pMinus=" << obj.pMinus << "]";
			return os;
		}

		~PpPm() {
        	delete pPlus;
        	delete pMinus;
    	}

	};

	static std::vector<PpPm> calcPpPm(MatrixCol matC) {
		std::vector<PpPm> result;
		for (int row = 0; row < matC.getRowCount(); row++) {
			result.push_back(PpPm(row));
		}
		for (int icol = 0, cole = matC.getColumnCount(); icol < cole; icol++) {
			SparseIntArray col = matC.getColumn(icol);
			for (int i = 0, ie = col.size(); i < ie; i++) {
				PpPm toedit = result[col.keyAt(i)];
				if (col.valueAt(i) < 0) {
					toedit.pMinus.append(icol, true);
				} else {
					toedit.pPlus.append(icol, true);
				}
			}
		}
		return result;
	}

	class Check11bResult {

	public:
		// The first columnindex where c_hj < 0 respectivly c_hj > 0
		const int col;
		// The whole row
		const int row;
		// The set P+ respectivly P-
		const SparseBoolArray p;

		/**
		 * Constructor to save the data.
		 *
		 * @paramk - the first column index where a component is less respectively
		 *         greater than zero.
		 * @param h     - the whole row.
		 * @param pPlus - the set P+ respectively P-.
		 */
		Check11bResult(int k, int row, SparseBoolArray pPlus) : col(k), row(row), pPlus(pPlus) {
		}
	};

//	static Check11bResult check11b(std::vector<PpPm> pppms, int startIndex) {
//		for (PpPm pppm : pppms.subList(startIndex, pppms.size())) {
//			Check11bResult res = check11bPppm(pppm);
//			if (res != nullptr) {
//				return res;
//			}
//		}
//		for (PpPm pppm : pppms.subList(0, startIndex)) {
//			Check11bResult res = check11bPppm(pppm);
//			if (res != null) {
//				return res;
//			}
//		}
//		return nullptr;
//	}
//
//	static Check11bResult check11bPppm(PpPm pppm) {
//		if (pppm.pMinus.size() == 1) {
//			return Check11bResult(pppm.pMinus.keyAt(0), pppm.row, pppm.pPlus);
//		} else if (pppm.pPlus.size() == 1) {
//			return Check11bResult(pppm.pPlus.keyAt(0), pppm.row, pppm.pMinus);
//		}
//		return nullptr;
//	}

}

#endif /* INVARIANTCALCULATOR_H_ */