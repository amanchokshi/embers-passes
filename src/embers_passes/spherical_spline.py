"""
spherical_spline.py
====================
JAX-native bivariate B-spline on R^2.

Coordinate convention
---------------------
This code is coordinate-agnostic: it is a plain tensor-product B-spline
over two real coordinates, here called ``u`` and ``v``. It has been used
both in native spherical coordinates (colatitude/azimuth) and in
orthographic projection coordinates on the unit disk without modification.
Nothing in the implementation assumes a spherical domain — any rectangular
(u, v) parameterization works.

Historical naming note: the class is called ``SphericalSpline`` and the
factory functions ``make_spherical_spline`` / ``make_spherical_spline_raw``
because this module's original motivating use case was a spline over the
upper hemisphere. The names stuck (and are used elsewhere in the codebase),
but what they build is a fully generic bivariate spline — the knot vector
fields are named ``t_u``/``t_v`` accordingly, not ``t_theta``/``t_phi``.

If you do use spherical coordinates directly, note that working right at
the pole can introduce a coordinate discontinuity; orthographic or another
locally well-behaved projection avoids this.

The spline is a tensor-product B-spline of degree (p, q) in (u, v).
All functions are pure JAX — jit-able, vmap-able, and differentiable
(grad / jacfwd / hessian), compatible with NumPyro / HMC.

Notes on differentiability
--------------------------
The Cox–de Boor recurrence produces piecewise-polynomial basis functions that
are p-1 times continuously differentiable in the interior of each knot span.
Differentiating through the degree-0 indicator step would introduce NaNs
from the Heaviside non-differentiability.  We therefore:

  * Use the **analytic B-spline derivative formula** (one degree lower basis)
    in ``eval_spline_gradient_3d`` / ``eval_spline_gradient_autodiff``.
  * Enable smooth jax.grad by replacing the evaluation function with a
    **de Boor algorithm** variant that avoids an explicit indicator: finding
    the knot span with ``jnp.searchsorted`` (non-differentiable w.r.t. x,
    but used only for index lookup, not in the arithmetic path).

The cleanest approach for NumPyro use: call ``eval_spline`` for forward
passes.  If you need ``jax.grad(eval_spline)`` directly, use the de Boor
evaluator variant ``eval_spline_deboor`` provided below.

Public API
----------
SphericalSpline               – pytree dataclass (generic bivariate spline)
make_spherical_spline         – from scipy RectBivariateSpline
make_spherical_spline_raw     – from numpy arrays directly
eval_spline                   – f(u,v) via basis-vector dot product
eval_spline_batch             – vmap over (N,) arrays
eval_spline_deboor            – f via de Boor (smooth jax.grad w.r.t. x)
eval_spline_deboor_batch      – vmap de Boor version
eval_spline_samples           – vmap eval_spline_batch over a leading
                                 sample axis of coefficients (S, n_u, n_v)
make_disk_mask                – boolean (n_u, n_v) mask of basis functions
                                 whose support intersects the unit disk
make_flat_index               – mask -> (flat_to_ij, n_active) index mapping
scatter_coeffs                – scatter n_active sampled coeffs into full
                                 (n_u, n_v) grid using a disk mask
make_constraint_info          – precompute pivot index/weights to enforce
                                 f(u0,v0) = f0 as a linear constraint
inject_pivot                  – insert the deterministic pivot coefficient
                                 into free coefficients given constraint info
make_knots_from_grid          – knot vectors from a rectangular (u,v) grid,
                                 via scipy's knot-placement heuristic
make_uniform_knots            – clamped uniform knot vectors, no scipy
                                 dependency
"""

from __future__ import annotations
import dataclasses
import jax
import jax.numpy as jnp
import numpy as np


# ---------------------------------------------------------------------------
# Register SphericalSpline as a JAX pytree so jit/vmap work naturally
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class SphericalSpline:
    """Tensor-product B-spline over two real coordinates (u, v).

    Named ``SphericalSpline`` for historical reasons (the original use case
    was a spline over the upper hemisphere), but this is a fully generic
    bivariate B-spline: u/v can be colatitude/azimuth, orthographic x/y,
    or any other rectangular parameterization.

    Attributes
    ----------
    t_u : (m_u,) knot vector in the first coordinate (u)
    t_v : (m_v,) knot vector in the second coordinate (v)
    c   : (n_u, n_v) B-spline coefficient matrix
    p   : degree in u (Python int, static)
    q   : degree in v (Python int, static)
    """
    t_u: jnp.ndarray
    t_v: jnp.ndarray
    c:   jnp.ndarray
    p:   int
    q:   int


# Register as pytree: arrays are leaves, (p,q) are auxiliary (static)
jax.tree_util.register_pytree_node(
    SphericalSpline,
    lambda spl: ([spl.t_u, spl.t_v, spl.c], (spl.p, spl.q)),
    lambda aux, children: SphericalSpline(*children, *aux),
)


# ---------------------------------------------------------------------------
# Factory helpers
# ---------------------------------------------------------------------------

def make_spherical_spline(rbs, dtype=None) -> SphericalSpline:
    """Convert ``scipy.interpolate.RectBivariateSpline`` → SphericalSpline.

    Parameters
    ----------
    rbs   : RectBivariateSpline fitted on a rectangular (u, v) grid.
    dtype : JAX dtype.  Defaults to jnp.float32; use jnp.float64 after
            ``jax.config.update("jax_enable_x64", True)``.
    """
    if dtype is None:
        dtype = jnp.float32
    tu  = np.array(rbs.tck[0])
    tv  = np.array(rbs.tck[1])
    ku  = int(rbs.degrees[0])
    kv  = int(rbs.degrees[1])
    nu  = len(tu) - ku - 1
    nv  = len(tv) - kv - 1
    c2d = np.array(rbs.tck[2]).reshape(nu, nv)
    return SphericalSpline(
        t_u = jnp.array(tu,  dtype=dtype),
        t_v   = jnp.array(tv,  dtype=dtype),
        c       = jnp.array(c2d, dtype=dtype),
        p=ku, q=kv,
    )


def make_spherical_spline_raw(t_u, t_v, c, p=3, q=3, dtype=None) -> SphericalSpline:
    """Build from numpy arrays (full knot vectors, coefficient matrix, degrees)."""
    if dtype is None:
        dtype = jnp.float32
    return SphericalSpline(
        t_u = jnp.array(t_u, dtype=dtype),
        t_v   = jnp.array(t_v, dtype=dtype),
        c       = jnp.array(c,   dtype=dtype),
        p=p, q=q,
    )


# ---------------------------------------------------------------------------
# B-spline basis via Cox–de Boor (static shapes, no dynamic slicing)
# ---------------------------------------------------------------------------
# Index arithmetic:
#   At step k, output has n_out = m-k-1 values (degree k, i=0..n_out-1).
#   Slices into the m-length knot vector:
#     t[i]     -> t[:n_out]
#     t[i+k]   -> t[k : k+n_out]
#     t[i+1]   -> t[1 : 1+n_out]
#     t[i+k+1] -> t[k+1 : k+1+n_out]   max index = k+1+n_out-1 = m-1  ✓

def _bspline_basis_active(t: jnp.ndarray, p: int, x) -> tuple:
    """The p+1 nonzero B-spline basis values at x, and the active span index.

    At any point x exactly p+1 basis functions B_{i*-p,p} ... B_{i*,p} are
    nonzero, where i* is the knot span satisfying t[i*] <= x < t[i*+1].
    This runs the triangular Cox–de Boor recurrence on a (p+1,) working array
    instead of the full (n,) basis vector, reducing work from O(n*p) to O(p^2).

    Parameters
    ----------
    t : (m,) knot vector with static shape
    p : degree (Python int, static at trace time)
    x : scalar (JAX-traced)

    Returns
    -------
    b      : (p+1,) nonzero basis values in order,
             b[j] = B_{i*-p+j, p}(x)  for j = 0 .. p
    i_star : scalar int — active span index, clamped to [p, n-1]
    """
    m  = t.shape[0]   # static
    n  = m - p - 1    # number of basis functions
    dt = t.dtype

    # Find active span: clamp to [p, n-1] so the window i*-p..i* is always valid.
    i_star = jnp.clip(jnp.searchsorted(t, x, side='right') - 1, p, n - 1)

    # Pull the knot window t[i*-p .. i*+p+1], length 2p+2 (static).
    # At recurrence step k, index j (0..p) accesses:
    #   ti   = t_win[j]       (t[i*-p+j])       max j   = p     < 2p+2 ✓
    #   tikl = t_win[j+k]     (t[i*-p+j+k])     max j+k = 2p    < 2p+2 ✓
    #   ti1  = t_win[j+1]     (t[i*-p+j+1])     max j+1 = p+1   < 2p+2 ✓
    #   tikr = t_win[j+k+1]   (t[i*-p+j+k+1])   max j+k+1=2p+1  < 2p+2 ✓
    t_win = jax.lax.dynamic_slice(t, (i_star - p,), (2 * p + 2,))

    # Initialise: degree-0 restricted to active span.
    # x lies in span i*, which is the last (index p) of our p+1 spans.
    b = jnp.zeros(p + 1, dtype=dt).at[p].set(1.0)

    js = jnp.arange(p + 1)   # static

    # Cox–de Boor recurrence, k = 1 .. p (Python loop → static unroll)
    for k in range(1, p + 1):
        ti   = t_win[js]           # t[i*-p+j]
        tikl = t_win[js + k]       # t[i*-p+j+k]
        ti1  = t_win[js + 1]       # t[i*-p+j+1]
        tikr = t_win[js + k + 1]   # t[i*-p+j+k+1]

        dl = tikl - ti
        dr = tikr - ti1
        al = jnp.where(dl != 0, (x - ti)   / dl, jnp.zeros(p + 1, dt))
        ar = jnp.where(dr != 0, (tikr - x) / dr, jnp.zeros(p + 1, dt))

        b_shifted = jnp.concatenate([b[1:], jnp.zeros(1, dtype=dt)])
        b = al * b + ar * b_shifted

    return b, i_star


# ---------------------------------------------------------------------------
# De Boor algorithm — smooth w.r.t. x for jax.grad
# ---------------------------------------------------------------------------
# The indicator in degree-0 basis is not differentiable.  The de Boor
# algorithm avoids this: the knot-span index i* = searchsorted(t, x) - 1 is
# found via a non-differentiable lookup (integer output), and then only smooth
# arithmetic follows.  jax.grad ignores the integer lookup and correctly
# differentiates the arithmetic.

def _deboor_eval(t: jnp.ndarray, p: int, c_1d: jnp.ndarray, x) -> jnp.ndarray:
    """Evaluate a 1-D B-spline via de Boor's algorithm (scalar, smooth in x).

    Parameters
    ----------
    t    : (m,) knot vector
    p    : degree (Python int)
    c_1d : (m-p-1,) coefficient vector
    x    : scalar evaluation point

    Returns
    -------
    scalar spline value
    """
    m = t.shape[0]
    n = m - p - 1   # number of basis functions = number of coefficients

    # Find knot span: i* such that t[i*] <= x < t[i*+1]
    # Use clipping to stay in [p, n-1] (valid interior spans)
    i_star = jnp.clip(jnp.searchsorted(t, x, side='right') - 1, p, n - 1)

    # Extract p+1 relevant coefficients (static window via dynamic_slice)
    # dynamic_slice requires static sizes; window size = p+1 is static.
    d = jax.lax.dynamic_slice(c_1d, (i_star - p,), (p + 1,))

    # De Boor triangular scheme (p rounds, each shrinking by 1)
    for r in range(1, p + 1):
        # j runs from i_star down to i_star - p + r
        # alpha[j] = (x - t[j]) / (t[j+p-r+1] - t[j])
        j_start = i_star - p + r  # j values: j_start .. i_star
        n_live  = p - r + 1       # number of d values still live
        js      = j_start + jnp.arange(n_live)
        t_j     = t[js]
        t_jp    = t[js + p - r + 1]
        denom   = t_jp - t_j
        alpha   = jnp.where(denom != 0, (x - t_j) / denom, 0.0)
        d = d.at[:n_live].set(
            (1.0 - alpha) * d[:n_live] + alpha * d[1 : n_live + 1]
        )

    return d[0]


def eval_spline_deboor(spl: SphericalSpline, u, v) -> jnp.ndarray:
    """Evaluate the spline via the de Boor algorithm (smooth → jax.grad works).

    Slightly more expensive than ``eval_spline`` but fully differentiable
    w.r.t. both u and v everywhere in the interior.

    Parameters
    ----------
    spl : SphericalSpline
    u   : scalar, first coordinate, within the domain of spl.t_u
    v   : scalar, second coordinate, within the domain of spl.t_v
    """
    # Evaluate in v for each u basis coefficient row, then in u
    n_u = spl.c.shape[0]
    n_v = spl.c.shape[1]

    # 1. For each u-row i, collapse the v dimension → v_vals[i]
    v_vals = jax.vmap(
        lambda c_row: _deboor_eval(spl.t_v, spl.q, c_row, v)
    )(spl.c)  # (n_u,)

    # 2. Evaluate the resulting 1-D spline in u
    return _deboor_eval(spl.t_u, spl.p, v_vals, u)


def eval_spline_deboor_batch(
    spl: SphericalSpline, u: jnp.ndarray, v: jnp.ndarray
) -> jnp.ndarray:
    """De Boor evaluation at N points → (N,) array."""
    return jax.vmap(lambda uu, vv: eval_spline_deboor(spl, uu, vv))(u, v)


# ---------------------------------------------------------------------------
# Spline evaluation (basis-vector version — fast forward pass)
# ---------------------------------------------------------------------------

def eval_spline(spl: SphericalSpline, u, v) -> jnp.ndarray:
    """Evaluate f(u, v) using only the (p+1)×(q+1) active coefficients.

    At any point only p+1 basis functions in u and q+1 in v are nonzero,
    so the bilinear form reduces to a dot product over a (p+1)×(q+1) subblock
    of c extracted with dynamic_slice.

    Use for forward passes.  For jax.grad w.r.t. (u, v), use
    ``eval_spline_deboor``.
    """
    b_u, i_u = _bspline_basis_active(spl.t_u, spl.p, u)  # (p+1,), scalar
    b_v, i_v = _bspline_basis_active(spl.t_v,   spl.q, v)  # (q+1,), scalar

    # Extract the (p+1)×(q+1) active subblock of c.
    # Slice start: i*-p is the first nonzero basis function index.
    c_sub = jax.lax.dynamic_slice(spl.c,
                                  (i_u - spl.p, i_v - spl.q),
                                  (spl.p + 1,   spl.q + 1))
    return jnp.dot(b_u, c_sub @ b_v)


def eval_spline_batch(
    spl: SphericalSpline, u: jnp.ndarray, v: jnp.ndarray
) -> jnp.ndarray:
    """Evaluate at N points → (N,) array."""
    return jax.vmap(lambda uu, vv: eval_spline(spl, uu, vv))(u, v)

def eval_spline_samples(
    spl: SphericalSpline,   # c has shape (S, n_u, n_v)
    u:   jnp.ndarray,       # scalar or (N,)
    v:   jnp.ndarray,       # scalar or (N,)
) -> jnp.ndarray:           # (S,) or (S, N)
    """Evaluate over S posterior samples of c.

    vmap is over the leading sample axis of spl.c; u/v are treated
    as fixed (in_axes=None broadcasts them across samples).
    """
    return jax.vmap(
        lambda c: eval_spline_batch(
            SphericalSpline(spl.t_u, spl.t_v, c, spl.p, spl.q),
            u, v,
        )
    )(spl.c)  # vmap sees spl.c as (S, n_u, n_v) and maps over axis 0

def make_disk_mask(t_u, t_v, p, q):
    """Boolean mask of shape (n_u, n_v): True where B_i(u)*B_j(v) has
    support intersecting the unit disk. Computed once at setup time (numpy).

    Strategy: the support of B_{i,p} is (t_u[i], t_u[i+p+1]).
    B_i(u)*B_j(v) can be nonzero on the disk iff the rectangle
    (t_u[i], t_u[i+p+1]) x (t_v[j], t_v[j+p+1]) intersects the unit disk.
    We check this geometrically: a rectangle intersects the unit disk iff
    the closest point in the rectangle to the origin has norm <= 1.
    """
    t_u = np.array(t_u); t_v = np.array(t_v)
    n_u = len(t_u) - p - 1
    n_v = len(t_v) - q - 1
    mask = np.zeros((n_u, n_v), dtype=bool)
    for i in range(n_u):
        u_lo, u_hi = t_u[i], t_u[i + p + 1]
        for j in range(n_v):
            v_lo, v_hi = t_v[j], t_v[j + q + 1]
            # Closest point in [u_lo,u_hi]x[v_lo,v_hi] to origin
            u_near = np.clip(0.0, u_lo, u_hi)
            v_near = np.clip(0.0, v_lo, v_hi)
            if u_near**2 + v_near**2 <= 1.0:
                mask[i, j] = True
    return mask

def make_flat_index(mask):
    """Map from flat sample index -> (i, j) coefficient index.
    Returns (flat_to_ij, n_active) where flat_to_ij has shape (n_active, 2).
    """
    ij = np.argwhere(mask)          # (n_active, 2)
    return ij, len(ij)

def scatter_coeffs(c_flat, mask, ij):
    """Scatter n_active sampled coefficients into full (n_u, n_v) array.
    JAX-compatible: uses jnp.zeros + indexed update (no Python control flow).

    Parameters
    ----------
    c_flat : (n_active,)  — the sampled coefficients
    mask   : (n_u, n_v)   — boolean mask (static, not traced)
    ij     : (n_active, 2) — integer indices (static)

    Returns
    -------
    c_full : (n_u, n_v)   — zero outside disk, c_flat inside
    """
    n_u, n_v = mask.shape
    c_full = jnp.zeros((n_u, n_v), dtype=c_flat.dtype)
    return c_full.at[ij[:, 0], ij[:, 1]].set(c_flat)

def make_constraint_info(t_u, t_v, p, q, u0, v0, ij):
    """
    Compute everything needed to enforce f(u0, v0) = f0 as a linear constraint
    by expressing one coefficient as a deterministic function of the others.

    At (u0, v0), only (p+1)*(q+1) basis functions are nonzero. Their indices
    in the full (n_u, n_v) grid are the 'active' set. We pick one of them
    (the 'pivot') to solve for deterministically, and the rest are sampled freely.

    Parameters
    ----------
    t_u, t_v : knot vectors (numpy)
    p, q     : degrees
    u0, v0   : reference point
    ij       : (n_active, 2) int array from make_flat_index — the mask-active indices
    f0       : desired function value at (u0, v0)

    Returns
    -------
    pivot_flat_idx : int
        Index into c_flat (the n_active-length vector) of the pivot coefficient.
        c_flat will have length n_active - 1 (free params); the pivot is inserted
        at this position.
    pivot_weight : float
        The basis value B_{i*,p}(u0) * B_{j*,q}(v0) for the pivot coefficient.
        The pivot coefficient = (f0 - sum of others * their weights) / pivot_weight.
    contrib_flat_idxs : (p+1)*(q+1) int array
        Indices into c_flat of ALL basis functions active at (u0,v0), including pivot.
    contrib_weights : (p+1)*(q+1) float array
        Corresponding basis values at (u0, v0).
    """
    t_u = np.array(t_u); t_v = np.array(t_v)
    n_u = len(t_u) - p - 1
    n_v = len(t_v) - q - 1

    # Active span indices at (u0, v0)
    i_star = int(np.clip(np.searchsorted(t_u, u0, side='right') - 1, p, n_u - 1))
    j_star = int(np.clip(np.searchsorted(t_v, v0, side='right') - 1, q, n_v - 1))

    # Global (i, j) indices of the (p+1)*(q+1) nonzero basis functions
    i_idxs = np.arange(i_star - p, i_star + 1)  # (p+1,)
    j_idxs = np.arange(j_star - q, j_star + 1)  # (q+1,)
    ii, jj = np.meshgrid(i_idxs, j_idxs, indexing='ij')
    ij_contrib = np.stack([ii.ravel(), jj.ravel()], axis=1)  # ((p+1)*(q+1), 2)

    # Evaluate basis values at (u0, v0) via de Boor (numpy, static)
    from scipy.interpolate import BSpline
    b_u = BSpline.design_matrix(np.array([u0]), t_u, p).toarray()[0]  # (n_u,)
    b_v = BSpline.design_matrix(np.array([v0]), t_v, q).toarray()[0]  # (n_v,)
    # Outer product, then extract the active block
    weights_2d = np.outer(b_u[i_idxs], b_v[j_idxs])  # (p+1, q+1)
    weights = weights_2d.ravel()                        # ((p+1)*(q+1),)

    # Map (i,j) -> flat index in c_flat via ij lookup
    # ij is (n_active, 2); build a dict for O(1) lookup
    ij_to_flat = {(r[0], r[1]): k for k, r in enumerate(ij)}
    contrib_flat_idxs = np.array([ij_to_flat[(r[0], r[1])] for r in ij_contrib])

    # Choose pivot: the contributor with the largest |weight| for numerical stability
    pivot_local = int(np.argmax(np.abs(weights)))
    pivot_flat_idx = int(contrib_flat_idxs[pivot_local])
    pivot_weight   = float(weights[pivot_local])

    return pivot_flat_idx, pivot_weight, contrib_flat_idxs, weights


def inject_pivot(c_free, pivot_flat_idx, pivot_weight,
                 contrib_flat_idxs, contrib_weights, f0):
    """
    Given n_active-1 free coefficients, insert the deterministic pivot so that
    f(u0, v0) = f0 is satisfied, returning a c_flat of length n_active.

    The constraint is:
        sum_k  contrib_weights[k] * c_flat[contrib_flat_idxs[k]]  =  f0
    Solving for the pivot:
        c_pivot = (f0 - sum_{k != pivot} w_k * c_k) / w_pivot

    Parameters
    ----------
    c_free           : (n_active - 1,) JAX array of free coefficients
    pivot_flat_idx   : int  (static)
    pivot_weight     : float (static)
    contrib_flat_idxs: (n_contrib,) int array (static)
    contrib_weights  : (n_contrib,) float array (static)
    f0               : desired value (scalar, may be JAX-traced)

    Returns
    -------
    c_flat : (n_active,) JAX array with pivot inserted
    """
    # Shift: c_free is indexed 0..n_active-2; positions >= pivot_flat_idx
    # in the final c_flat are at position-1 in c_free.
    # Insert a placeholder zero at pivot position, then overwrite.
    c_flat = jnp.concatenate([
        c_free[:pivot_flat_idx],
        jnp.zeros(1, dtype=c_free.dtype),
        c_free[pivot_flat_idx:],
    ])

    # Indices of non-pivot contributors in c_flat (after insertion, so same
    # as contrib_flat_idxs since pivot slot was zero before overwrite)
    mask_not_pivot = contrib_flat_idxs != pivot_flat_idx
    other_idxs    = contrib_flat_idxs[mask_not_pivot]   # static
    other_weights = contrib_weights[mask_not_pivot]      # static

    # Compute pivot value
    other_sum  = jnp.dot(jnp.array(other_weights), c_flat[other_idxs])
    c_pivot    = (f0 - other_sum) / pivot_weight

    return c_flat.at[pivot_flat_idx].set(c_pivot)


def make_knots_from_grid(
    u_grid: np.ndarray,
    v_grid: np.ndarray,
    p:      int = 3,
    q:      int = 3,
    dtype=None,
):
    """Pre-compute B-spline knot vectors from a rectangular grid.

    Call this ONCE outside any JAX-traced function.  The returned knot
    vectors are fixed for the lifetime of the model; only the coefficient
    matrix ``c`` needs to be traced / inferred.

    Uses scipy's ``RectBivariateSpline`` internals to place knots optimally
    for the given grid, but discards the fit itself.  Alternatively you can
    supply uniform knots directly.

    Parameters
    ----------
    u_grid : (n_u,) 1-D grid of first-coordinate values
    v_grid : (n_v,) 1-D grid of second-coordinate values
    p, q   : B-spline degrees (default 3 = cubic)
    dtype  : JAX dtype for the returned arrays (default jnp.float32;
             use jnp.float64 after enabling x64)

    Returns
    -------
    t_u : jnp.ndarray, shape (n_u + p + 1,)  – knot vector in u
    t_v : jnp.ndarray, shape (n_v + q + 1,)  – knot vector in v
    p   : int
    q   : int
    n_u : int  – number of basis functions / coefficients in u
    n_v : int  – number of basis functions / coefficients in v

    Example
    -------
    >>> t_u, t_v, p, q, n_u, n_v = make_knots_from_grid(
    ...     u_grid, v_grid, p=3, q=3, dtype=jnp.float64)
    >>> # Inside NumPyro model:
    >>> c   = numpyro.sample("c", dist.Normal(0, 1).expand([n_u, n_v]))
    >>> spl = SphericalSpline(t_u, t_v, c, p, q)
    >>> val = eval_spline_deboor(spl, u_obs, v_obs)
    """
    from scipy.interpolate import RectBivariateSpline

    if dtype is None:
        dtype = jnp.float32

    # Fit a dummy spline just to get scipy's knot placement
    dummy_z = np.zeros((len(u_grid), len(v_grid)))
    rbs     = RectBivariateSpline(u_grid, v_grid, dummy_z, kx=p, ky=q)

    t_u_np = np.array(rbs.tck[0])
    t_v_np = np.array(rbs.tck[1])
    n_u    = len(t_u_np) - p - 1
    n_v    = len(t_v_np) - q - 1

    t_u = jnp.array(t_u_np, dtype=dtype)
    t_v = jnp.array(t_v_np, dtype=dtype)

    return t_u, t_v, p, q, n_u, n_v


def make_uniform_knots(
    u_min: float,
    u_max: float,
    v_min: float,
    v_max: float,
    n_u:   int,
    n_v:   int,
    p:     int = 3,
    q:     int = 3,
    dtype=None,
):
    """Build clamped uniform knot vectors (no scipy dependency).

    Produces knot vectors with ``p+1`` repeated boundary knots (clamped
    B-spline) and ``n_u - p - 1`` interior knots uniformly spaced.

    Parameters
    ----------
    u_min/max : domain bounds in the first coordinate
    v_min/max : domain bounds in the second coordinate
    n_u       : number of basis functions in u
    n_v       : number of basis functions in v
    p, q      : degrees

    Returns
    -------
    Same as ``make_knots_from_grid``.
    """
    if dtype is None:
        dtype = jnp.float32

    def clamped_knots(a, b, n, deg):
        n_interior = n - deg - 1
        interior   = np.linspace(a, b, n_interior + 2)[1:-1]  # exclude endpoints
        return np.concatenate([
            np.repeat(a, deg + 1),
            interior,
            np.repeat(b, deg + 1),
        ])

    t_u_np = clamped_knots(u_min, u_max, n_u, p)
    t_v_np = clamped_knots(v_min, v_max, n_v, q)

    t_u = jnp.array(t_u_np, dtype=dtype)
    t_v = jnp.array(t_v_np, dtype=dtype)

    return t_u, t_v, p, q, n_u, n_v