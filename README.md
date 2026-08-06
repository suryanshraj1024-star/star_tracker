# Star Tracker — Algorithm README

A four-phase pipeline that ingests raw 2D sensor coordinates, maps them into 3D celestial space, identifies stars against a reference catalog using angular invariance, and then solves for spacecraft attitude via SVD-based rotation estimation with noise-robustness validation.

---

## Core Programmatic Logic Paths

### Phase 1 — 2D Sensor → 3D Unit Vectors

Raw pixel/sensor coordinates `(x, y)` are lifted into 3D by treating the focal length `f` as the z-component, then L2-normalized to produce a unit direction vector on the celestial sphere.

```
normalize(x, y, f):
    mag = sqrt(x² + y² + f²)
    return (x/mag, y/mag, f/mag)
```

Each observed star becomes a unit vector `v̂ ∈ ℝ³` encoding only direction, not distance.

### Phase 2 — Angular Separation (Rotation-Invariant Feature)

For every unique pair of observed stars, the inter-star angle is computed via the dot product:

```
angular_sep(a, b):
    dot = clip(a · b, −1, 1)
    return (dot, arccos(dot))
```

The same computation is performed exhaustively over all catalog pairs. Because angular separations are invariant to spacecraft rotation, they serve as the matching fingerprint between observed sky and reference catalog.

### Phase 3 — Star Identification via Voting

For each observed pair whose angular separation is within tolerance `tol = 0.005 rad` of a catalog pair, both catalog stars accumulate a vote for both observed stars. After all pair comparisons are done, a greedy one-to-one assignment is resolved: each observed star is assigned to its highest-voted, as-yet-unassigned catalog star. This ensures a valid bijective mapping with no catalog star used twice.

### Phase 4 — Attitude Estimation (SVD Rotation Solver) and Noise Robustness

Using the matched observed-to-catalog vector pairs, the spacecraft's rotation matrix `R` is solved via Singular Value Decomposition on the cross-covariance matrix `H`:

```
H = Σ obs_i ⊗ cat_i    (outer product sum)
U, S, Vt = SVD(H)
R = U @ Vt
```

A determinant check enforces a proper rotation (`det(R) = +1`), correcting reflections if they arise. The resulting `R` is decomposed into ZYX Euler angles (yaw/pitch/roll). Noise robustness is then evaluated by injecting Gaussian pixel noise at escalating sigma levels (0.05 mm → 10.0 mm) over 30 Monte Carlo trials each, with a pass threshold of ≥ 80% correct identification rate.

---

## Spatial Data Mapping Architecture

```
Sensor Plane (2D)          Celestial Sphere (3D)         Catalog Frame (3D)
  (x, y) px coords    →    unit vector v̂ = (vx, vy, vz)  ←  pre-computed catalog vectors
      ↓                              ↓                              ↓
  focal length f           inter-star angles (Δθ)          catalog pair angles
  as z-coordinate          rotation-invariant feature      same invariant feature
                                     ↓                              ↓
                              Angular Matching (|Δθ_obs − Δθ_cat| ≤ tol)
                                          ↓
                               Vote-weighted bijective assignment
                                          ↓
                           Cross-covariance matrix H ∈ ℝ³ˣ³
                                          ↓
                              SVD → Rotation Matrix R
                                          ↓
                              ZYX Euler Angles (yaw, pitch, roll)
```

Two coordinate frames flow through the pipeline:

- **Sensor/Camera frame** — pixel or millimeter coordinates on a flat detector, origin at the optical center.
- **Celestial reference frame** — unit vectors on the unit sphere, origin at the spacecraft/camera origin. The catalog lives here in a fixed inertial frame; the observed vectors live here in the current body frame. `R` is the rotation that maps one frame to the other.

The focal plane visualization (Phase 4's 3D plot) makes this explicit: observed LOS rays project through the focal plane at `z = f_scene`, and catalog stars are plotted on the unit sphere, with matched ones highlighted.

---

## Matrix Optimization Logic — Analytical Explanation

The core optimization problem is: *given a set of matched vector pairs `{(b̂_i, r̂_i)}`* (body-frame observed vs. reference-frame catalog), find the rotation matrix `R ∈ SO(3)` that minimizes the total alignment error:

```
minimize  Σᵢ ‖R · b̂ᵢ − r̂ᵢ‖²
subject to  Rᵀ R = I,  det(R) = +1
```

This is Wahba's Problem — a classical constrained least-squares problem on the rotation group. The closed-form solution proceeds in three steps:

**Step 1 — Cross-covariance matrix.** Construct `H = B Rᵀ` in aggregate as:

```
H = Σᵢ b̂ᵢ rᵀᵢ   ∈ ℝ³ˣ³
```

`H` captures the directional correlation between the two point sets.

**Step 2 — SVD factorization.** Decompose `H = U Σ Vᵀ`. The singular values in `Σ` measure how well the two frames align along each principal axis; `U` and `V` encode the respective orientations of those axes in body and reference frame.

**Step 3 — Optimal rotation.** The rotation that minimizes Wahba's loss is:

```
R* = U diag(1, 1, det(U)·det(V)) Vᵀ
```

The diagonal correction term (equivalent to flipping the sign of the last column of `Vt` when `det < 0`) enforces `det(R) = +1`, projecting the SVD solution from the full orthogonal group `O(3)` onto the proper rotation group `SO(3)`. Without it, the optimizer can find a reflection (an improper rotation with determinant −1) that yields a lower algebraic loss but is physically meaningless as a spacecraft attitude.

The approach is optimal in the least-squares sense, closed-form (no iterative solver needed), numerically stable for small star counts, and directly produces Euler angles for downstream attitude control.

---

## Dependencies

| Library | Role |
|---|---|
| `numpy` | Vector math, SVD, rotation algebra |
| `matplotlib` | 3D celestial sphere visualization |
| `mpl_toolkits.mplot3d` | 3D axes and polygon collections |
| `random` / `collections` | Monte Carlo noise trials, vote tallying |
| `itertools` | Combinatorial star-pair enumeration |
