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


class InvariantCalculator {

	static const boolean DEBUG = false;

	/**
	 * Hidden constructor
	 */
	private InvariantCalculator() {
	}

	/**
	 * Enumeration for choosing which algorithm should be used.
	 */
	public enum InvariantAlgorithm {

		FARKAS, PIPE;
	}

	/**
	 * A class for holding the sets P+ = {j | c_hj &gt; 0} and P- = {j | c_hj less
	 * 0} for a given row.
	 */
	private static class PpPm {

		// The row
		public final int row;
		// P+ set
		public final SparseBoolArray pPlus = new SparseBoolArray();
		// P- set
		public final SparseBoolArray pMinus = new SparseBoolArray();

		/**
		 * initially empty.
		 *
		 * @param row
		 */
		public PpPm(int row) {
			this.row = row;
		}

		public void setValue(int j, int val) {
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

		@Override
		public String toString() {
			return "PpPm [row=" + row + ", pPlus=" + pPlus + ", pMinus=" + pMinus + "]";
		}

	}

	/**
	 * Calculates for a given matrix the P+ = {j | c_hj &gt; 0} and P- = {j | c_hj
	 * less 0} sets for each row.
	 *
	 * @param matC - the matrix from which the sets should be calculated.
	 * @return The result of the calculation
	 */
	private static List<PpPm> calcPpPm(IntMatrixCol matC) {
		final List<PpPm> result = new ArrayList<>();
		for (int row = 0; row < matC.getRowCount(); row++) {
			result.add(new PpPm(row));
		}
		for (int icol = 0, cole = matC.getColumnCount(); icol < cole; icol++) {
			SparseIntArray col = matC.getColumn(icol);
			for (int i = 0, ie = col.size(); i < ie; i++) {
				PpPm toedit = result.get(col.keyAt(i));
				if (col.valueAt(i) < 0) {
					toedit.pMinus.append(icol, true);
				} else {
					toedit.pPlus.append(icol, true);
				}
			}
		}
		return result;
	}

	/**
	 * Holds the result of the check11b-check. That means it holds the row, the
	 * column index where a component is less than respectively greater than zero
	 * and the P+ respectively P- set if there exists a row in the given matrix such
	 * that |P+| == 1 or |P-| == 1.
	 */
	private static class Check11bResult {

		// The first columnindex where c_hj < 0 respectivly c_hj > 0
		public final int col;
		// The whole row
		public final int row;
		// The set P+ respectivly P-
		public final SparseBoolArray p;

		/**
		 * Constructor to save the data.
		 *
		 * @paramk - the first column index where a component is less respectively
		 *         greater than zero.
		 * @param h     - the whole row.
		 * @param pPlus - the set P+ respectively P-.
		 */
		public Check11bResult(int k, int row, SparseBoolArray pPlus) {
			this.col = k;
			this.row = row;
			this.p = pPlus;
		}
	}

	/**
	 * Checks if there exists a row in the given matrix such that |P+| == 1 or |P-|
	 * == 1 and returns the row or null if such a row do not exists.
	 *
	 * @param pppms      the list of all rows with P+ and P- sets.
	 * @param startIndex
	 * @return the row which satisfy |P+| == 1 or |P-| == 1 or null if not existent.
	 */
	private static Check11bResult check11b(List<PpPm> pppms, int startIndex) {
		for (PpPm pppm : pppms.subList(startIndex, pppms.size())) {
			Check11bResult res = check11bPppm(pppm);
			if (res != null) {
				return res;
			}
		}
		for (PpPm pppm : pppms.subList(0, startIndex)) {
			Check11bResult res = check11bPppm(pppm);
			if (res != null) {
				return res;
			}
		}
		return null;
	}

	private static Check11bResult check11bPppm(PpPm pppm) {
		if (pppm.pMinus.size() == 1) {
			return new Check11bResult(pppm.pMinus.keyAt(0), pppm.row, pppm.pPlus);
		} else if (pppm.pPlus.size() == 1) {
			return new Check11bResult(pppm.pPlus.keyAt(0), pppm.row, pppm.pMinus);
		}
		return null;
	}

	/**
	 * Calculates the invariants with the algorithm based on
	 * http://pipe2.sourceforge.net/documents/PIPE-Report.pdf (page 19).
	 *
	 * @param mat          - the matrix to calculate the invariants from.
	 * @param onlyPositive whether we just stop at Flows or go for Semi-Flows
	 * @param pnames       variable names
	 * @return a generator set of the invariants.
	 */
	public static Set<SparseIntArray> calcInvariantsPIPE(IntMatrixCol mat, boolean onlyPositive) {
		if (mat.getColumnCount() == 0 || mat.getRowCount() == 0) {
			return new HashSet<>();
		}
		IntMatrixCol tmat = mat.transpose();
		Set<SparseIntArray> normed = new HashSet<>();
		for (int i = 0; i < tmat.getColumnCount(); i++) {
			SparseIntArray norm = tmat.getColumn(i);
			normalize(norm);
			normed.add(norm);
		}
		if (normed.size() < tmat.getColumnCount()) {
			System.out.println("Normalized transition count is " + normed.size() + " out of " + tmat.getColumnCount()
					+ " initially.");
		}
		IntMatrixCol matnorm = new IntMatrixCol(tmat.getRowCount(), 0);
		for (SparseIntArray col : normed) {
			matnorm.appendColumn(col);
		}
		final IntMatrixCol matB = phase1PIPE(matnorm.transpose());

//		final MatrixCol matB = phase1PIPE(new MatrixCol(mat));
		// We want to work with columns in this part of the algorithm
		// We add and remove columns all day => we want to switch to a column based
		// representation
		// order of rows is really irrelevant + columns which are identical up to
		// scaling factor are useless
		// let's use a set of columns.
		Set<SparseIntArray> colsBsparse = new HashSet<>(2 * matB.getColumnCount());
		for (int i = 0; i < matB.getColumnCount(); i++) {
			SparseIntArray col = matB.getColumn(i);
			if (col.size() != 0) {
				normalizeWithSign(col);
				colsBsparse.add(col);
			}
		}

		if (!onlyPositive) {
			return colsBsparse;
		}
	}

#endif /* INVARIANTCALCULATOR_H_ */