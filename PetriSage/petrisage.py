#!/usr/bin/env python3
import sys
import logging
from sage.all import Matrix, ZZ, QQ, lcm, vector
from sage.libs.pari import pari

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger('PetriSage')

def read_sparse_matrix(filename, transpose=False, as_pari=False):
    logger.info(f"Reading {filename}, transpose={transpose}, as_pari={as_pari}")
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    rows, cols = map(int, lines[0].split())
    logger.info(f"Header: {rows} rows, {cols} cols")
    entries = {}
    for line in lines[1:]:
        parts = line.split()
        if not parts:
            continue
        row = int(parts[0])
        for pair in parts[1:]:
            col, val = map(int, pair.split(':'))
            r, c = (col, row) if transpose else (row, col)
            if r >= (cols if transpose else rows) or c >= (rows if transpose else cols):
                logger.warning(f"Out-of-range entry: ({r},{c}) for {rows}x{cols}")
                continue
            entries[(r, c)] = val
    dims = (cols, rows, entries) if transpose else (rows, cols, entries)
    logger.info(f"Loaded matrix: {dims[0]}x{dims[1]}, {len(entries)} non-zero entries")
    logger.debug(f"Entries: {entries}")
    if as_pari:
        m, n = dims[0], dims[1]
        # Create a dense PARI matrix with 0-based indexing.
        pari_mat = pari.matrix(m, n)
        for (r, c), val in entries.items():
            pari_mat[r, c] = val  # 0-based indexing for cypari2
        logger.debug(f"PARI matrix:\n{pari_mat}")
        return m, n, pari_mat
    return dims

def compute_flows_hnf(m, n, sparse_data):
    logger.info(f"Computing HNF and right kernel for {m}x{n} matrix")
    C = Matrix(ZZ, m, n, sparse=True, entries=sparse_data)
    logger.debug(f"Matrix C:\n{C}")
    H, U = C.hermite_form(transformation=True)
    logger.debug(f"HNF H:\n{H}")
    logger.debug(f"Transformation U:\n{U}")
    rank = C.rank()
    nullity = n - rank
    logger.info(f"Rank: {rank}, Nullity: {nullity}")
    K = C.right_kernel()
    flows = K.basis()
    logger.debug(f"Right kernel basis:\n{flows}")
    logger.info(f"Extracted {len(flows)} flows")
    logger.debug(f"Flows: {flows}")
    return flows

def compute_flows_pari_kernel(m, n, pari_mat):
    logger.info(f"Computing PARI kernel for {m}x{n} matrix")
    logger.debug(f"PARI matrix:\n{pari_mat}")
    # Compute the kernel via PARI's global function.
    K_pari = pari.matker(pari_mat)
    logger.debug(f"PARI kernel matrix (raw): {K_pari}")
    # Convert the PARI matrix into a Sage matrix.
    # Note: We iterate over rows; pari.matker returns a matrix in row form.
    K_pari_list = []
    for row in list(K_pari):
        K_pari_list.append([int(x) for x in list(row)])
    sage_K = Matrix(ZZ, K_pari_list)
    logger.debug(f"Sage kernel matrix (before transpose):\n{sage_K}")
    # Transpose so that nonzero kernel basis vectors become columns.
    sage_K = sage_K.transpose()
    logger.debug(f"Sage kernel matrix (after transpose):\n{sage_K}")
    # Filter out any columns that are identically zero.
    flows = []
    for i in range(sage_K.ncols()):
        col = sage_K.column(i)
        if any(coef != 0 for coef in col):
            flows.append(col)
    logger.info(f"Extracted {len(flows)} flows")
    logger.debug(f"Flows: {flows}")
    return flows

def compute_flows_snf(m, n, sparse_data):
    logger.info(f"Computing SNF for {m}x{n} matrix")
    C = Matrix(ZZ, m, n, sparse=True, entries=sparse_data)
    logger.debug(f"Matrix C:\n{C}")
    D, S, T = C.smith_form()
    logger.debug(f"SNF D:\n{D}")
    logger.debug(f"Left transformation S:\n{S}")
    logger.debug(f"Right transformation T:\n{T}")
    
    # Compute rank from D (number of non-zero diagonal entries)
    rank = sum(1 for i in range(min(m, n)) if D[i, i] != 0)
    nullity = n - rank
    logger.info(f"Rank: {rank}, Nullity: {nullity}")
    
    # Extract kernel basis from the last (n - rank) columns of T
    if nullity > 0:
        flows = [T.column(i) for i in range(n - nullity, n)]
    else:
        flows = []
    
    logger.info(f"Extracted {len(flows)} flows directly from SNF")
    logger.debug(f"Flows: {flows}")
    
    # Optional: Verify flows are in the kernel (for debugging)
    for flow in flows:
        assert C * flow == vector(ZZ, m), "Flow not in kernel!"
    
    return flows

def compute_flows_snf_old(m, n, sparse_data):
    logger.info(f"Computing SNF for {m}x{n} matrix")
    C = Matrix(ZZ, m, n, sparse=True, entries=sparse_data)
    logger.debug(f"Matrix C:\n{C}")
    D, _, T = C.smith_form()
    logger.debug(f"SNF D:\n{D}")
    logger.debug(f"Transformation T:\n{T}")
    rank = sum(1 for i in range(min(m, n)) if D[i, i] != 0)
    nullity = n - rank
    logger.info(f"Rank: {rank}, Nullity: {nullity}")
    K = C.right_kernel()
    flows = K.basis()
    logger.debug(f"Right kernel basis:\n{flows}")
    logger.info(f"Extracted {len(flows)} flows")
    logger.debug(f"Flows: {flows}")
    return flows

def compute_flows_rational(m, n, sparse_data):
    logger.info(f"Computing rational kernel for {m}x{n} matrix")
    C = Matrix(QQ, m, n, sparse=True, entries=sparse_data)
    logger.debug(f"Matrix C (QQ):\n{C}")
    K = C.right_kernel()
    basis = K.basis()
    flows = []
    for v in basis:
        denom = lcm([x.denominator() for x in v if x != 0])
        flows.append((denom * v).change_ring(ZZ))
    logger.debug(f"Rational basis:\n{basis}")
    logger.info(f"Extracted {len(flows)} flows")
    logger.debug(f"Flows: {flows}")
    return flows

def write_flows(flows, outfile, mode='TFLOWS'):
    logger.info(f"Writing {len(flows)} flows to {outfile}")
    with open(outfile, 'w') as f:
        f.write(f"{len(flows)}\n")
        for flow in flows:
            terms = [(coeff, i + 1) for i, coeff in enumerate(flow) if coeff != 0]
            f.write(f"{len(terms)}")
            for coeff, idx in terms:
                f.write(f" {coeff} {idx}")
            f.write("\n")
        f.write("0\n")

def print_usage():
    print("""
Usage: petrisage.py <model.mat> <output_file> <TFLOWS|PFLOWS> [--backend=BACKEND]
Computes P- or T-flows (integer invariants) of a Petri net incidence matrix using SageMath.

Arguments:
  <model.mat>     Input file with sparse incidence matrix.
  <output_file>   Output file for flows in GreatSPN .tba/.pba format.
  <TFLOWS|PFLOWS> Mode: TFLOWS for transition flows, PFLOWS for place flows.

Options:
  --backend=BACKEND  Select computation backend (default: hnf):
    hnf            Hermite Normal Form via FLINT (sparse, integer basis).
    pari_kernel    Direct integer kernel via PARI/GP (dense matrix).
    snf            Smith Normal Form via FLINT/PARI (integer basis).
    rational       Rational kernel via LinBox, scaled to integers.
""")

def main():
    if len(sys.argv) < 4 or sys.argv[1] in ['-h', '--help']:
        print_usage()
        sys.exit(1 if len(sys.argv) < 4 else 0)
    matrix_file = sys.argv[1]
    output_file = sys.argv[2]
    mode = sys.argv[3].upper()
    backend = 'hnf'
    for arg in sys.argv[4:]:
        if arg.startswith('--backend='):
            backend = arg[len('--backend='):].lower()
    if mode not in ['TFLOWS', 'PFLOWS']:
        logger.error("Mode must be TFLOWS or PFLOWS")
        sys.exit(1)
    backends = {
        'hnf': compute_flows_hnf,
        'pari_kernel': compute_flows_pari_kernel,
        'snf': compute_flows_snf,
        'rational': compute_flows_rational
    }
    if backend not in backends:
        logger.error(f"Unknown backend '{backend}'. Available: {', '.join(backends.keys())}")
        sys.exit(1)
    as_pari = (backend == 'pari_kernel')
    m, n, data = read_sparse_matrix(matrix_file, transpose=(mode == 'PFLOWS'), as_pari=as_pari)
    flows = backends[backend](m, n, data)
    write_flows(flows, output_file, mode)
    logger.info(f"Computed {len(flows)} {mode.lower()[:-1]}s using {backend}, saved to {output_file}")

if __name__ == "__main__":
    main()
