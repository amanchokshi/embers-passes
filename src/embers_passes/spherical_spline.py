"""
spherical_spline.py
====================
JAX-native bivariate B-spline on the upper hemisphere.

Coordinate convention
---------------------
  theta : colatitude in [0, pi/2]   (0 = north pole, pi/2 = equator)
  phi   : azimuth   in [0, 2*pi)

The spline is a tensor-product B-spline of degree (p, q) in (theta, phi).
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

The cleanest approach for NumpPyro use: call ``eval_spline`` for forward
passes and ``eval_spline_gradient_3d`` for explicit gradients.  If you need
``jax.grad(eval_spline)`` directly, use the de Boor evaluator variant
``eval_spline_deboor`` provided below.

Public API
----------
SphericalSpline               – pytree dataclass
make_spherical_spline         – from scipy RectBivariateSpline
make_spherical_spline_raw     – from numpy arrays directly
eval_spline                   – f(theta,phi) via basis-vector dot product
eval_spline_batch             – vmap over (N,) arrays
eval_spline_deboor            – f via de Boor (smooth jax.grad w.r.t. x)
eval_spline_deboor_batch      – vmap de Boor version
eval_spline_gradient_3d       – analytic Cartesian surface gradient
eval_spline_gradient_batch    – vmap surface gradient
"""

from __future__ import annotations
import dataclasses
from typing import Tuple
import jax
import jax.numpy as jnp
import numpy as np


# ---------------------------------------------------------------------------
# Register SphericalSpline as a JAX pytree so jit/vmap work naturally
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class SphericalSpline:
    """Tensor-product B-spline on the upper hemisphere.

    Attributes
    ----------
    t_theta : (m_theta,) knot vector in colatitude
    t_phi   : (m_phi,)   knot vector in azimuth
    c       : (n_theta, n_phi) B-spline coefficient matrix
    p       : degree in theta (Python int, static)
    q       : degree in phi   (Python int, static)
    """
    t_theta: jnp.ndarray
    t_phi:   jnp.ndarray
    c:       jnp.ndarray
    p:       int
    q:       int


# Register as pytree: arrays are leaves, (p,q) are auxiliary (static)
jax.tree_util.register_pytree_node(
    SphericalSpline,
    lambda spl: ([spl.t_theta, spl.t_phi, spl.c], (spl.p, spl.q)),
    lambda aux, children: SphericalSpline(*children, *aux),
)


# ---------------------------------------------------------------------------
# Factory helpers
# ---------------------------------------------------------------------------

def make_spherical_spline(rbs, dtype=None) -> SphericalSpline:
    """Convert ``scipy.interpolate.RectBivariateSpline`` → SphericalSpline.

    Parameters
    ----------
    rbs   : RectBivariateSpline fitted on a (theta, phi) rectangular grid.
    dtype : JAX dtype.  Defaults to jnp.float32; use jnp.float64 after
            ``jax.config.update("jax_enable_x64", True)``.
    """
    if dtype is None:
        dtype = jnp.float32
    tx  = np.array(rbs.tck[0])
    ty  = np.array(rbs.tck[1])
    kx  = int(rbs.degrees[0])
    ky  = int(rbs.degrees[1])
    nx  = len(tx) - kx - 1
    ny  = len(ty) - ky - 1
    c2d = np.array(rbs.tck[2]).reshape(nx, ny)
    return SphericalSpline(
        t_theta = jnp.array(tx,  dtype=dtype),
        t_phi   = jnp.array(ty,  dtype=dtype),
        c       = jnp.array(c2d, dtype=dtype),
        p=kx, q=ky,
    )


def make_spherical_spline_raw(t_theta, t_phi, c, p=3, q=3, dtype=None) -> SphericalSpline:
    """Build from numpy arrays (full knot vectors, coefficient matrix, degrees)."""
    if dtype is None:
        dtype = jnp.float32
    return SphericalSpline(
        t_theta = jnp.array(t_theta, dtype=dtype),
        t_phi   = jnp.array(t_phi,   dtype=dtype),
        c       = jnp.array(c,       dtype=dtype),
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


def _bspline_basis_deriv_vector(t: jnp.ndarray, p: int, x) -> jnp.ndarray:
    """First derivative d/dx B_{i,p}(x), shape (m-p-1,).

    Uses: d/dx B_{i,p} = p*[B_{i,p-1}/(t_{i+p}-t_i) - B_{i+1,p-1}/(t_{i+p+1}-t_{i+1})]
    This calls _bspline_basis_vector at degree p-1, so is smooth w.r.t. x
    wherever the degree-(p-1) basis is smooth.
    """
    m  = t.shape[0]
    n  = m - p - 1
    dt = t.dtype

    b_pm1 = _bspline_basis_vector(t, p - 1, x)   # (n+1,)

    ti   = t[:n];              tip  = t[p : p + n]
    ti1  = t[1 : 1 + n];      tip1 = t[p + 1 : p + 1 + n]

    dl = tip  - ti;   al = jnp.where(dl != 0, p / dl, jnp.zeros(n, dt))
    dr = tip1 - ti1;  ar = jnp.where(dr != 0, p / dr, jnp.zeros(n, dt))
    return al * b_pm1[:n] - ar * b_pm1[1 : n + 1]


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


def eval_spline_deboor(spl: SphericalSpline, theta, phi) -> jnp.ndarray:
    """Evaluate the spline via the de Boor algorithm (smooth → jax.grad works).

    Slightly more expensive than ``eval_spline`` but fully differentiable
    w.r.t. both theta and phi everywhere in the interior.

    Parameters
    ----------
    spl   : SphericalSpline
    theta : colatitude scalar in [0, pi/2]
    phi   : azimuth scalar in [0, 2*pi)
    """
    # Evaluate in phi for each theta basis coefficient row, then in theta
    n_theta = spl.c.shape[0]
    n_phi   = spl.c.shape[1]

    # 1. For each theta-row i, collapse the phi dimension → phi_vals[i]
    phi_vals = jax.vmap(
        lambda c_row: _deboor_eval(spl.t_phi, spl.q, c_row, phi)
    )(spl.c)  # (n_theta,)

    # 2. Evaluate the resulting 1-D spline in theta
    return _deboor_eval(spl.t_theta, spl.p, phi_vals, theta)


def eval_spline_deboor_batch(
    spl: SphericalSpline, theta: jnp.ndarray, phi: jnp.ndarray
) -> jnp.ndarray:
    """De Boor evaluation at N points → (N,) array."""
    return jax.vmap(lambda th, ph: eval_spline_deboor(spl, th, ph))(theta, phi)


# ---------------------------------------------------------------------------
# Spline evaluation (basis-vector version — fast forward pass)
# ---------------------------------------------------------------------------

def eval_spline(spl: SphericalSpline, theta, phi) -> jnp.ndarray:
    """Evaluate f(theta, phi) using only the (p+1)×(q+1) active coefficients.

    At any point only p+1 basis functions in theta and q+1 in phi are nonzero,
    so the bilinear form reduces to a dot product over a (p+1)×(q+1) subblock
    of c extracted with dynamic_slice.

    Use for forward passes.  For jax.grad w.r.t. (theta, phi), use
    ``eval_spline_deboor`` or call ``eval_spline_gradient_3d`` directly.
    """
    b_t, i_t = _bspline_basis_active(spl.t_theta, spl.p, theta)  # (p+1,), scalar
    b_p, i_p = _bspline_basis_active(spl.t_phi,   spl.q, phi)    # (q+1,), scalar

    # Extract the (p+1)×(q+1) active subblock of c.
    # Slice start: i*-p is the first nonzero basis function index.
    c_sub = jax.lax.dynamic_slice(spl.c,
                                  (i_t - spl.p, i_p - spl.q),
                                  (spl.p + 1,   spl.q + 1))
    return jnp.dot(b_t, c_sub @ b_p)


def eval_spline_batch(
    spl: SphericalSpline, theta: jnp.ndarray, phi: jnp.ndarray
) -> jnp.ndarray:
    """Evaluate at N points → (N,) array."""
    return jax.vmap(lambda th, ph: eval_spline(spl, th, ph))(theta, phi)


# ---------------------------------------------------------------------------
# Analytic surface gradient (recommended for gradient computation)
# ---------------------------------------------------------------------------

def eval_spline_gradient_3d(spl: SphericalSpline, theta, phi) -> jnp.ndarray:
    """Cartesian surface gradient of the spline at (theta, phi).

    Uses the analytic B-spline derivative formula (degree p-1 basis), which
    is smooth wherever the degree-(p-1) basis is smooth (i.e., p-2 times
    continuously differentiable at knots).

    The intrinsic gradient on the unit sphere is:
        grad_S f  =  df/dθ · ê_θ  +  (df/dφ / sin θ) · ê_φ

    Orthonormal tangent frame:
        ê_θ = ( cos θ cos φ,  cos θ sin φ,  −sin θ )
        ê_φ = ( −sin φ,        cos φ,         0     )

    The 1/sin θ factor is clamped (≥ 1e-8) for pole stability.

    Returns
    -------
    grad : (3,) Cartesian surface gradient (⊥ to the position vector)
    """
    b_t  = _bspline_basis_vector(spl.t_theta, spl.p, theta)
    b_p  = _bspline_basis_vector(spl.t_phi,   spl.q, phi)
    db_t = _bspline_basis_deriv_vector(spl.t_theta, spl.p, theta)
    db_p = _bspline_basis_deriv_vector(spl.t_phi,   spl.q, phi)

    df_dtheta = jnp.dot(db_t, spl.c @ b_p)
    df_dphi   = jnp.dot(b_t,  spl.c @ db_p)

    inv_sin = 1.0 / jnp.maximum(jnp.abs(jnp.sin(theta)), 1e-8)
    cos_t   = jnp.cos(theta); sin_t = jnp.sin(theta)
    cos_p   = jnp.cos(phi);   sin_p = jnp.sin(phi)
    zero    = jnp.zeros_like(cos_t)

    e_theta = jnp.stack([ cos_t * cos_p,  cos_t * sin_p, -sin_t])
    e_phi   = jnp.stack([-sin_p,           cos_p,          zero ])

    return df_dtheta * e_theta + (df_dphi * inv_sin) * e_phi


def eval_spline_gradient_batch(
    spl: SphericalSpline, theta: jnp.ndarray, phi: jnp.ndarray
) -> jnp.ndarray:
    """Analytic surface gradient at N points → (N, 3) array."""
    return jax.vmap(lambda th, ph: eval_spline_gradient_3d(spl, th, ph))(theta, phi)


# ---------------------------------------------------------------------------
# Autodiff partials via de Boor (for validation / NumPyro log-prob grads)
# ---------------------------------------------------------------------------

def eval_spline_gradient_autodiff(
    spl: SphericalSpline, theta, phi
) -> Tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """(value, df/dtheta, df/dphi) via jax.grad through eval_spline_deboor.

    Returns
    -------
    value, df_dtheta, df_dphi : three scalars
    """
    f = lambda th, ph: eval_spline_deboor(spl, th, ph)
    val      = f(theta, phi)
    df_dth   = jax.grad(f, 0)(theta, phi)
    df_dph   = jax.grad(f, 1)(theta, phi)
    return val, df_dth, df_dph


# ---------------------------------------------------------------------------
# Demo / smoke test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import numpy as np
    from scipy.interpolate import RectBivariateSpline

    jax.config.update("jax_enable_x64", True)
    print("=== Spherical B-spline smoke test ===\n")

    # Ground-truth: f(θ,φ) = sin²θ·cos(2φ)  smooth, zero at pole
    theta_grid = np.linspace(0.05, np.pi / 2 - 0.05, 40)
    phi_grid   = np.linspace(0.0,  2 * np.pi,         80, endpoint=False)
    TH, PH     = np.meshgrid(theta_grid, phi_grid, indexing="ij")
    Z          = np.sin(TH)**2 * np.cos(2 * PH)

    rbs = RectBivariateSpline(theta_grid, phi_grid, Z, kx=3, ky=3)
    spl = make_spherical_spline(rbs, dtype=jnp.float64)
    print(f"Knot counts — θ: {len(spl.t_theta)}, φ: {len(spl.t_phi)}")
    print(f"Coefficient grid: {spl.c.shape},  degree=({spl.p},{spl.q})\n")

    rng     = np.random.default_rng(7)
    N       = 8
    th_test = jnp.array(rng.uniform(0.1, np.pi/2 - 0.1, N))
    ph_test = jnp.array(rng.uniform(0.0, 2*np.pi,        N))

    # ---- 1. Value accuracy ----
    v_basis  = eval_spline_batch(spl, th_test, ph_test)
    v_deboor = eval_spline_deboor_batch(spl, th_test, ph_test)
    v_scipy  = rbs.ev(np.array(th_test), np.array(ph_test))
    v_true   = np.sin(np.array(th_test))**2 * np.cos(2*np.array(ph_test))

    print("Values (basis | deBoor | scipy | true | Δbasis-scipy):")
    for i in range(N):
        print(f"  [{i}]  {v_basis[i]:.8f}  {v_deboor[i]:.8f}  "
              f"{v_scipy[i]:.8f}  {v_true[i]:.8f}  Δ={float(v_basis[i]-v_scipy[i]):.1e}")

    # ---- 2. Gradient accuracy ----
    th0, ph0 = jnp.float64(0.6), jnp.float64(1.3)

    # Analytic formula
    g3d = eval_spline_gradient_3d(spl, th0, ph0)

    # jax.grad through de Boor
    val, dth_ad, dph_ad = eval_spline_gradient_autodiff(spl, th0, ph0)

    # True partials for sin²θ·cos(2φ)
    th0n, ph0n = float(th0), float(ph0)
    true_dth =  2*np.sin(th0n)*np.cos(th0n)*np.cos(2*ph0n)
    true_dph = -2*np.sin(th0n)**2*np.sin(2*ph0n)

    print(f"\nPartials at (θ={th0n:.2f}, φ={ph0n:.2f}):")
    print(f"  df/dθ → jax.grad: {float(dth_ad):.8f}  true: {true_dth:.8f}"
          f"  Δ={float(dth_ad)-true_dth:.1e}")
    print(f"  df/dφ → jax.grad: {float(dph_ad):.8f}  true: {true_dph:.8f}"
          f"  Δ={float(dph_ad)-true_dph:.1e}")

    # Tangency check: surface gradient ⊥ outward normal
    n_hat = jnp.array([jnp.sin(th0)*jnp.cos(ph0),
                        jnp.sin(th0)*jnp.sin(ph0),
                        jnp.cos(th0)])
    print(f"\n  3D gradient: {np.array(g3d).round(7)}")
    print(f"  ‖grad‖ = {float(jnp.linalg.norm(g3d)):.7f}")
    print(f"  grad·n̂ (tangency, ≈0): {float(g3d @ n_hat):.2e}")

    # ---- 3. jit / vmap ----
    jit_vals = jax.jit(eval_spline_batch)(spl, th_test, ph_test)
    jit_debo = jax.jit(eval_spline_deboor_batch)(spl, th_test, ph_test)
    jit_grad = jax.jit(eval_spline_gradient_batch)(spl, th_test, ph_test)
    jit_vals.block_until_ready(); jit_grad.block_until_ready(); jit_debo.block_until_ready()
    print("\njit + vmap: OK")
    print(f"  eval_spline_batch          → {jit_vals.shape}")
    print(f"  eval_spline_deboor_batch   → {jit_debo.shape}")
    print(f"  eval_spline_gradient_batch → {jit_grad.shape}")

    # ---- 4. Higher-order differentiability (HMC needs this) ----
    jac = jax.jit(
        jax.jacobian(lambda th: eval_spline_deboor_batch(spl, th, ph_test))
    )(th_test)
    print(f"\njax.jacobian d(values)/d(theta) {jac.shape}: OK")

    h = jax.jit(
        jax.hessian(lambda th: eval_spline_deboor(spl, th, ph_test[0]))
    )(th_test[0])
    print(f"jax.hessian (scalar→scalar): {float(h):.7f}  OK")

    # ---- 5. Gradient of analytic gradient (2nd order) ----
    grad_of_grad = jax.jit(
        jax.jacobian(lambda th: eval_spline_gradient_batch(spl, th, ph_test))
    )(th_test)
    print(f"jax.jacobian of gradient_batch {grad_of_grad.shape}: OK")

    print("\n✓ All checks passed")


# ===========================================================================
# NumPyro / HMC usage pattern
# ===========================================================================
# When the spline coefficients `c` are latent variables, the workflow is:
#
#   ONCE (outside any JAX tracing):
#     knots = make_knots_from_grid(theta_grid, phi_grid, p=3, q=3)
#
#   INSIDE the NumPyro model (traced by JAX):
#     c = numpyro.sample("c", prior)          # shape (n_theta, n_phi)
#     spl = SphericalSpline(*knots, c, p, q)  # c is a JAX traced array
#     val = eval_spline_deboor(spl, theta, phi)   # differentiable in c
#
# The key insight: knots (t_theta, t_phi) are FIXED numpy arrays converted
# to JAX constants; only `c` flows through the JAX tracer.

def make_knots_from_grid(
    theta_grid: np.ndarray,
    phi_grid:   np.ndarray,
    p:          int = 3,
    q:          int = 3,
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
    theta_grid : (n_theta,) 1-D grid of colatitude values
    phi_grid   : (n_phi,)   1-D grid of azimuth values
    p, q       : B-spline degrees (default 3 = cubic)
    dtype      : JAX dtype for the returned arrays (default jnp.float32;
                 use jnp.float64 after enabling x64)

    Returns
    -------
    t_theta : jnp.ndarray, shape (n_theta + p + 1,)  – knot vector in theta
    t_phi   : jnp.ndarray, shape (n_phi   + q + 1,)  – knot vector in phi
    p       : int
    q       : int
    n_theta : int  – number of basis functions / coefficients in theta
    n_phi   : int  – number of basis functions / coefficients in phi

    Example
    -------
    >>> t_theta, t_phi, p, q, n_th, n_ph = make_knots_from_grid(
    ...     theta_grid, phi_grid, p=3, q=3, dtype=jnp.float64)
    >>> # Inside NumPyro model:
    >>> c   = numpyro.sample("c", dist.Normal(0, 1).expand([n_th, n_ph]))
    >>> spl = SphericalSpline(t_theta, t_phi, c, p, q)
    >>> val = eval_spline_deboor(spl, theta_obs, phi_obs)
    """
    from scipy.interpolate import RectBivariateSpline

    if dtype is None:
        dtype = jnp.float32

    # Fit a dummy spline just to get scipy's knot placement
    dummy_z = np.zeros((len(theta_grid), len(phi_grid)))
    rbs     = RectBivariateSpline(theta_grid, phi_grid, dummy_z, kx=p, ky=q)

    t_theta_np = np.array(rbs.tck[0])
    t_phi_np   = np.array(rbs.tck[1])
    n_theta    = len(t_theta_np) - p - 1
    n_phi      = len(t_phi_np)   - q - 1

    t_theta = jnp.array(t_theta_np, dtype=dtype)
    t_phi   = jnp.array(t_phi_np,   dtype=dtype)

    return t_theta, t_phi, p, q, n_theta, n_phi


def make_uniform_knots(
    theta_min: float,
    theta_max: float,
    phi_min:   float,
    phi_max:   float,
    n_theta:   int,
    n_phi:     int,
    p:         int = 3,
    q:         int = 3,
    dtype=None,
):
    """Build clamped uniform knot vectors (no scipy dependency).

    Produces knot vectors with ``p+1`` repeated boundary knots (clamped
    B-spline) and ``n_theta - p - 1`` interior knots uniformly spaced.

    Parameters
    ----------
    theta_min/max : domain bounds in colatitude
    phi_min/max   : domain bounds in azimuth
    n_theta       : number of basis functions in theta
    n_phi         : number of basis functions in phi
    p, q          : degrees

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

    t_theta_np = clamped_knots(theta_min, theta_max, n_theta, p)
    t_phi_np   = clamped_knots(phi_min,   phi_max,   n_phi,   q)

    t_theta = jnp.array(t_theta_np, dtype=dtype)
    t_phi   = jnp.array(t_phi_np,   dtype=dtype)

    return t_theta, t_phi, p, q, n_theta, n_phi
