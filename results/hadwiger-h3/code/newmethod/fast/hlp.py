"""Thin direct HiGHS call that reproduces scipy.optimize.linprog(method='highs') on
min c.x s.t. A_ub x <= b_ub, lb <= x <= ub  -- same matrix layout (CSC of the dense A), same options
(presolve on, dual simplex strategy, no output), so results are bit-identical, without linprog's
per-call input cleaning / option checking / sparse construction / dual extraction overhead."""
import numpy as np
from scipy.optimize._highspy import _core as _h
from scipy.optimize._highspy._core import simplex_constants as _sc
INF = _h.kHighsInf
_OPT = _h.HighsOptions()
_OPT.presolve = 'on'
_OPT.highs_debug_level = _h.HighsDebugLevel.kHighsDebugLevelNone
_OPT.log_to_console = False
_OPT.output_flag = False
_OPT.simplex_strategy = _sc.SimplexStrategy.kSimplexStrategyDual
_OPTIMAL = _h.HighsModelStatus.kOptimal

def csc_dense(A):
    """CSC arrays of a dense matrix, entries in column-major order, zeros dropped (as scipy csc_array(A))."""
    At = A.T
    nz = At != 0
    cnt = nz.sum(axis=1)
    indptr = np.zeros(A.shape[1] + 1, dtype=np.int32); np.cumsum(cnt, out=indptr[1:])
    cols, rows = np.nonzero(nz)
    return indptr, rows.astype(np.int32), At[cols, rows]

def solve_ub(c, A, b, lb, ub, csc=None):
    """returns (ok, x) ; ok iff HiGHS model status is optimal."""
    m, n = A.shape if csc is None else (len(b), len(c))
    indptr, indices, data = csc_dense(A) if csc is None else csc
    lp = _h.HighsLp()
    lp.num_col_ = n; lp.num_row_ = m
    lp.a_matrix_.num_col_ = n; lp.a_matrix_.num_row_ = m
    lp.a_matrix_.format_ = _h.MatrixFormat.kColwise
    lp.col_cost_ = c; lp.col_lower_ = lb; lp.col_upper_ = ub
    lp.row_lower_ = np.full(m, -INF); lp.row_upper_ = b
    lp.a_matrix_.start_ = indptr; lp.a_matrix_.index_ = indices; lp.a_matrix_.value_ = data
    hs = _h._Highs()
    hs.passOptions(_OPT)
    if hs.passModel(lp) == _h.HighsStatus.kError: return False, None
    if hs.run() == _h.HighsStatus.kError: return False, None
    if hs.getModelStatus() != _OPTIMAL: return False, None
    return True, np.array(hs.getSolution().col_value)

def bounds_arrays(bounds):
    lb = np.array([-INF if l is None else float(l) for l, u in bounds])
    ub = np.array([INF if u is None else float(u) for l, u in bounds])
    return lb, ub
