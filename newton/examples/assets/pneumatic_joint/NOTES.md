# Pneumatic joint example — asset notes & de-risking results

Scaffolded ahead of implementing `example_softbody_pneumatic_joint.py` (Wenlong/Gaozhang 2024
variable-stiffness soft joint). All findings below were **verified by running Newton 1.14.0 / Warp
1.14.0 locally** (RTX 3060, CUDA 12.9) — not assumptions.

## Assets present
- `chamber_right.msh`, `chamber_left.msh` — gmsh tet meshes of the two semicircular chambers
  (from `project files/.../Objects Felix/`). Units: **mm**. Right chamber: 1055 tets, 357 verts,
  extent **68.6 × 28 × 75 mm** (semicircular cross-section in X–Y, chamber axis along **Z**).
- `joint_bending_ground_truth.csv` — 3624 rows: `Time, Joint_angle[deg], chamber_1_pressure_bar,
  chamber_2_pressure_bar`. Verified: angle 0–39.86°, P1 0–0.742 bar, P2 0–0.355 bar, no NaNs,
  corr(angle,P1)=0.80. **bar → Pa = ×1e5.**
- `material_params.json`, `joint_params.py` — material fit + paper parameters.
- ⏳ **MISSING — `hinge.stl`**: Felix to export `IJ_Vertical_Hinge*` from SolidWorks.

## Verified API path (Phase A — works end to end)
```python
tm = newton.TetMesh.create_from_file(".../chamber_right.msh")
tm.custom_attributes = {}        # REQUIRED: gmsh files carry a 'gmsh:physical' custom attribute;
                                 # add_soft_mesh raises unless it is cleared or registered first.
builder.add_soft_mesh(pos=..., rot=..., scale=0.001, vel=wp.vec3(0.0),   # scale 0.001 = mm -> m
                      mesh=tm, density=1070, k_mu=3.97e5, k_lambda=..., k_damp=1.0,
                      validate_mesh=True)
builder.color()                  # REQUIRED before finalize() for SolverVBD
model = builder.finalize()
solver = newton.solvers.SolverVBD(model=model, iterations=...)   # VBD step confirmed to run
```
Combined 2-chamber model: **713 particles, 2109 tets, 1410 surface tris**. A single VBD step runs on
CPU and CUDA.

### TetMesh attribute names (the docstring's "indices" is wrong)
`vertices` (N,3) · `tet_indices` (4·tets, flat) · `surface_tri_indices` (3·tris, flat) ·
`tet_count` · `vertex_count` · `density/k_mu/k_lambda/k_damp` (all **None** for these meshes —
must be supplied) · `custom_attributes`.

## Material (from material_params.json)
- `k_mu` = mu_full = **3.97e5 Pa** (matches paper NH μ̄≈0.207–0.4 MPa).
- `k_lambda` (full, ν≈0.5): K − ⅔μ ≈ **2.08e7 Pa** (near-incompressible). VBD ran one step fine with
  this, but **sustained pressure stability is still unproven** — keep a softened λ (effective ν≈0.45)
  as a fallback knob. ρ = 1070 kg/m³.

## Phase B — cavity-face classification (the pressure-kernel input)
Pressure must act on the **interior cavity** faces only. Measured on the right chamber (706 surface
tris): inner and outer **arch walls separate cleanly by radius** (inner ≈ 31 mm, outer ≈ 35 mm from the
arch axis at (cx,cy)=(34.29, 21.0)). **A single radius threshold is NOT enough** — the flat chord faces
(Y≈28 outer / Y≈32 inner) and the 4 mm end-caps span many radii and pollute the inner set.
Robust classifier = radius band around R_in **AND** inward-pointing normal, **plus** the inner-chord
plane and interior-cap faces. Alternative: reuse Felix's existing gmsh boundary groups in
`project files/.../*_groups_bnd.pos`. This is real (but tractable) Phase-B work.

## Phase B — pressure kernel PROTOTYPED & validated
Prototype: `project files/prototype_cavity_pressure.py` (single free chamber, runs on CUDA).
- **Cavity classification via connected components works perfectly:** the enclosed void is a
  separate surface component. Right chamber → exactly 2 components: exterior shell (418 tris,
  max_r=35.0=R_out) and **cavity (288 tris, max_r=31.0=R_in)**. No thresholds needed.
- **Follower pressure kernel works:** `apply_cavity_pressure` adds `P · area_vec / 3` per cavity
  triangle into `particle_f` each substep (recomputed from live positions = true follower load).
  **Bake triangle winding outward ONCE at rest** — do NOT flip per-step against a fixed centroid
  (that injects energy as the body deforms → blew up to ±600 mm in testing). Direction confirmed
  correct: pins on bottom cap, top cap extends +Z.
- **CRITICAL FINDING — fibre reinforcement is mandatory.** A bare unreinforced silicone shell
  **balloons radially and inverts/NaNs** even at 0.03 bar, no matter the damping. This matches the
  paper's premise (and the buckling discussion in §4.3). Approximating the reinforcement with surface
  **membrane stiffness** (`tri_ke=tri_ka=5e4`, `tri_kd=1e2`) + heavy `k_damp=1e4` + velocity clamp
  (`particle_max_velocity=1.0`) gave a **stable, inversion-free run for all 60 frames** at 0.05 bar
  — but underdamped (oscillatory). **Next:** a proper reinforcement model (e.g. stiff circumferential
  hoop/fibre constraints) + quasi-static settling (slower ramp, more substeps, hold-to-settle) are
  needed for a physical bending result. This is the #1 modeling task before validation.

## Still open (unchanged from plan risks)
- Soft↔rigid end-cap coupling to the hinge (no native particle→body attach).
- Fibre reinforcement (axial-only extension) not represented by the bare mesh.
- Chosen meshes are R=35/L=75 mm, not the paper's R=19/L=49 mm → MVP validation is qualitative.
