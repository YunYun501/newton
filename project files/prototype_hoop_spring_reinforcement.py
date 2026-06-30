# SPDX-License-Identifier: Apache-2.0
"""
Pneumatic variable-stiffness joint (Gaozhang 2024) — antagonistic bending.

Variant: HOOP-SPRING fibre reinforcement. The inextensible fibre layer of the real
chamber is modelled by stiff CIRCUMFERENTIAL (hoop) springs added between existing wall
nodes of each chamber mesh — 1-D spring elements ringing the cross-section, NOT a
separate body or shell. They are stiff in the hoop direction (so the section cannot
balloon out or buckle/collapse inward) but leave the axial direction free (so pressure
still extends the chamber lengthwise -> bending). SolverVBD evaluates the springs
implicitly, so unlike an explicit penalty force they do not destabilise the solve.

TWO SEPARATE semicircular (D-shaped) silicone FEM chambers, one in each slot either
side of the central rigid hinge wall (they do NOT join across the middle). The single
+x D-mesh is instanced as-is for the right slot and rotated 180 deg about the length
axis for the left slot. Inflating one chamber extends that side axially, so the joint
bends about the pin; inflating the other bends it back. Co-inflating both stiffens
with little net bend (the antagonistic principle).

  * follower cavity-pressure load (custom Warp kernel) per cavity
  * bottom caps pinned (hard projection pin); top free -> the joint bends
  * quasi-static mass scaling + no ground plane keep the tiny-tet VBD solve stable
  * bending angle theta read out from top-cap vs base-cap centroids

Run (opens an OpenGL window):
  uv run --extra examples python "project files/prototype_cavity_pressure.py"
"""

import json
import os
from collections import defaultdict

import numpy as np
import warp as wp

import newton
import newton.examples

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CADDIR = os.path.join(ROOT, "project files", "CAD_files")
MESH = os.path.join(CADDIR, "chamber_half_thick.msh")
MATJSON = os.path.join(ROOT, "project files", "data", "material_params.json")

# Hinge geometry (mm). The two STLs are each exported about their own origin, so
# hinge_top must be raised by HINGE_TOP_DZ so the two pin knuckles coincide (forming
# the hinge) and the flanges land ~90 mm apart, matching the simplified-joint assembly.
# The revolute/bending axis is then +y through (x=0, z=PIN_Z_MM).
HINGE_BOTTOM = os.path.join(CADDIR, "hinge_bottom.stl")
HINGE_TOP = os.path.join(CADDIR, "hinge_top.stl")
HINGE_TOP_DZ_MM = 33.8
PIN_Z_MM = 34.6


GROUND_TRUTH = os.path.join(ROOT, "project files", "data", "ground_truth_experimental_data.csv")


def load_mesh(path):
    import trimesh  # noqa: PLC0415
    m = trimesh.load(path)
    return newton.Mesh(np.asarray(m.vertices, dtype=np.float32),
                       np.asarray(m.faces, dtype=np.int32).reshape(-1))


def load_ground_truth(path):
    """Experiment-1 log: joint angle [deg] + chamber 1/2 pressures [bar] over time."""
    d = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    return (np.radians(d["Joint_angle"].astype(float)),     # rad
            d["chamber_1_pressure_bar"].astype(float) * 1e5,  # Pa
            d["chamber_2_pressure_bar"].astype(float) * 1e5)  # Pa


# ---------------------------------------------------------------------------
# Geometry auto-detection + cavity classification.
# ---------------------------------------------------------------------------
def analyze_geometry(v):
    ext = v.max(0) - v.min(0)
    axis = int(np.argmax(ext))
    others = [i for i in range(3) if i != axis]
    cx, cy = v[:, others[0]].mean(), v[:, others[1]].mean()
    return {
        "axis": axis, "others": others, "center": (cx, cy),
        "L": float(ext[axis]), "z_min": float(v[:, axis].min()), "z_max": float(v[:, axis].max()),
        "r_out": float(np.hypot(v[:, others[0]] - cx, v[:, others[1]] - cy).max()),
        "scale": 0.001 if ext[axis] > 1.0 else 1.0,  # >1 native unit => millimetres
    }


def _surface_components(tris):
    n = len(tris)
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    edge_map = {}
    for t, (a, b, c) in enumerate(tris):
        for e in ((a, b), (b, c), (a, c)):
            key = (min(e), max(e))
            if key in edge_map:
                ra, rb = find(edge_map[key]), find(t)
                if ra != rb:
                    parent[ra] = rb
            else:
                edge_map[key] = t
    comps = defaultdict(list)
    for t in range(n):
        comps[find(t)].append(t)
    return list(comps.values())


def classify_cavities(v, surface_tris, geom):
    """Return a list of (tris, centroid) for each interior cavity surface.

    A fully-enclosed cavity is its own connected component of the boundary mesh;
    the outer shell is the component reaching the largest radius. Everything else
    is a cavity (this body has two).
    """
    tris = np.asarray(surface_tris).reshape(-1, 3)
    groups = _surface_components(tris)
    o0, o1 = geom["others"]
    cx, cy = geom["center"]
    rad = np.hypot(v[tris][:, :, o0].mean(1) - cx, v[tris][:, :, o1].mean(1) - cy)
    outer = int(np.argmax([rad[g].max() for g in groups]))
    cavities = []
    for i, g in enumerate(groups):
        if i == outer:
            continue
        ct = tris[np.array(g)]
        cavities.append((ct, v[np.unique(ct)].mean(0)))
    print(f"  {len(groups)} surface components -> {len(cavities)} cavities "
          f"(outer shell {len(groups[outer])} tris dropped)")
    return cavities


def bake_outward(tris, v, centroid):
    """Reorder each triangle so its winding normal points away from the cavity
    centroid (= the inflate direction), once at rest."""
    out = tris.copy()
    for t, (a, b, c) in enumerate(tris):
        nrm = np.cross(v[b] - v[a], v[c] - v[a])
        if np.dot(nrm, v[[a, b, c]].mean(0) - centroid) < 0.0:
            out[t] = [a, c, b]
    return out


def signed_cavity_volume(tris, q):
    """Enclosed volume of a closed, outward-wound triangle surface (divergence theorem):
    V = 1/6 * sum( a . (b x c) ). With ``tris`` baked outward, V>0 = a healthy cavity;
    V shrinking under pressure = the wall is being driven INWARD (implosion)."""
    a = q[tris[:, 0]]
    b = q[tris[:, 1]]
    c = q[tris[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def build_hoop_springs(builder, node_idx, pos, axis_xy, band_h, ke, kd, max_len):
    """Add inextensible-fibre hoop springs along a chamber's wall cross-section.

    The wall surface nodes are binned into thin bands along the length axis. Within each
    band the inner-cavity ring and the outer-wall ring (split by radius) are ordered by
    angle and ADJACENT nodes are tied by a short spring. The cross-section is a half-D
    (an open arc + a flat chord), so the rings are NOT wrapped closed and any pair farther
    apart than ``max_len`` is skipped — this keeps every spring local to the wall and
    avoids long springs jumping across the cavity (which would distort the mesh). The
    open hoops still resist the arc lengthening (ballooning); the axial direction is free
    so pressure still extends the chamber -> bending.

    Args:
        node_idx: global particle indices of the chamber's wall (surface) nodes.
        pos: their world rest positions, shape [n, 3].
        axis_xy: (cx, cy) of the chamber axis in the cross-section plane.
        band_h: band thickness along the length axis [m].
        max_len: skip any spring longer than this [m] (drops cross-cavity jumps).
    """
    z = pos[:, 2]
    rad = np.hypot(pos[:, 0] - axis_xy[0], pos[:, 1] - axis_xy[1])
    ang = np.arctan2(pos[:, 1] - axis_xy[1], pos[:, 0] - axis_xy[0])
    nb = max(1, int(np.ceil((z.max() - z.min()) / band_h)))
    edges = np.linspace(z.min(), z.max() + 1e-9, nb + 1)
    count = 0

    def link_arc(members):
        nonlocal count
        if len(members) < 3:
            return
        order = members[np.argsort(ang[members])]
        for k in range(len(order) - 1):       # open arc, no wrap-around
            i, j = order[k], order[k + 1]
            if np.linalg.norm(pos[i] - pos[j]) < max_len:
                builder.add_spring(int(node_idx[i]), int(node_idx[j]), ke, kd, 0.0)
                count += 1

    for b in range(nb):
        band = np.where((z >= edges[b]) & (z < edges[b + 1]))[0]
        if len(band) < 4:
            continue
        mid = np.median(rad[band])
        link_arc(band[rad[band] < mid])   # inner cavity arc
        link_arc(band[rad[band] >= mid])  # outer wall arc
    return count


# ---------------------------------------------------------------------------
# Warp kernels
# ---------------------------------------------------------------------------
@wp.kernel
def apply_cavity_pressure(
    particle_q: wp.array[wp.vec3],
    cavity_tris: wp.array[wp.vec3i],
    pressure: wp.array[float],
    particle_f: wp.array[wp.vec3],
):
    tid = wp.tid()
    tri = cavity_tris[tid]
    i, j, k = tri[0], tri[1], tri[2]
    area_vec = 0.5 * wp.cross(particle_q[j] - particle_q[i], particle_q[k] - particle_q[i])
    f = pressure[0] * area_vec / 3.0
    wp.atomic_add(particle_f, i, f)
    wp.atomic_add(particle_f, j, f)
    wp.atomic_add(particle_f, k, f)


@wp.kernel
def apply_fibre(
    particle_q: wp.array[wp.vec3],
    above: wp.array[int],       # 1 if the node sits above the pin (rest)
    rest_rad: wp.array[float],  # rest distance from the chamber centreline [m]
    pin_z: float,
    theta: wp.array[float],
    k_fibre: float,
    particle_f: wp.array[wp.vec3],
):
    # Inextensible fibre layer: hold each node at its REST radial distance from the
    # chamber centreline (free axially and tangentially) so cavity pressure cannot balloon
    # the wall radially and is forced into AXIAL extension. The centreline is vertical
    # below the pin and tilts by the hinge angle above it, so a rigid bend preserves the
    # radius (no resistance to bending) while local ballooning is opposed.
    tid = wp.tid()
    p = particle_q[tid]
    c = wp.cos(theta[0])
    s = wp.sin(theta[0])
    if above[tid] == 1:
        # un-rotate about the pin by -theta (inverse of the project_top R_y(-theta))
        rx = p[0]
        rz = p[2] - pin_z
        ux = rx * c + rz * s
        uy = p[1]
        r = wp.sqrt(ux * ux + uy * uy)
        if r > 1.0e-6:
            fmag = -k_fibre * (r - rest_rad[tid])
            fux = fmag * ux / r
            fuy = fmag * uy / r
            # rotate the radial force back into the bent frame (R_y(-theta))
            wp.atomic_add(particle_f, tid, wp.vec3(fux * c, fuy, fux * s))
    else:
        r = wp.sqrt(p[0] * p[0] + p[1] * p[1])
        if r > 1.0e-6:
            fmag = -k_fibre * (r - rest_rad[tid])
            wp.atomic_add(particle_f, tid, wp.vec3(fmag * p[0] / r, fmag * p[1] / r, 0.0))


@wp.kernel
def pin_particles(
    pin_idx: wp.array[int],
    pin_rest: wp.array[wp.vec3],
    relax: float,
    vdamp: float,
    particle_q: wp.array[wp.vec3],
    particle_qd: wp.array[wp.vec3],
):
    # Compliant projection pin. A HARD snap-to-rest (relax=1, vdamp=0) every substep makes
    # the boundary infinitely stiff: the first FREE tet layer above the bond starts each
    # solve already sheared and inverts (the "black hole"). Under-relaxing the correction
    # (relax<1) leaves pinned nodes partway, so the boundary tet is never fully sheared,
    # while accumulating over the substeps it still anchors the base. vdamp keeps some
    # velocity instead of zeroing it, removing the per-substep velocity discontinuity.
    tid = wp.tid()
    p = pin_idx[tid]
    particle_q[p] = particle_q[p] + relax * (pin_rest[tid] - particle_q[p])
    particle_qd[p] = vdamp * particle_qd[p]


# --- Emergent floating-rigid-plate hinge (pure FEM, angle from node displacement) ---
# Both chamber tops are silicone-bonded to ONE rigid distal plate that pivots about the
# central hinge. We model it by shape-matching each substep: fit the best rigid transform
# (rotation theta about +y PLUS an x-z translation) of the combined top-cap nodes, then
# project them onto it. The translation lets the plate centroid follow the chambers as they
# EXTEND -- a fixed-distance-to-pin rotation would freeze the chamber length and buckle the
# wall. Rotating about an OFF-CENTRE hinge pin is exactly (rotate about centroid) + (translate
# the centroid along the pin arc), so the x-z translation IS the off-centre-pivot behaviour.
# The bend angle theta EMERGES from the differential extension -- nothing is prescribed.
# Out-of-plane motion (y translation, x/z rotation) stays locked -> a clean single-axis hinge.
@wp.kernel
def accumulate_centroid(particle_q: wp.array[wp.vec3], top_idx: wp.array[int],
                        csum: wp.array[float]):
    tid = wp.tid()
    p = top_idx[tid]
    wp.atomic_add(csum, 0, particle_q[p][0])
    wp.atomic_add(csum, 1, particle_q[p][1])
    wp.atomic_add(csum, 2, particle_q[p][2])


@wp.kernel
def accumulate_fit(
    particle_q: wp.array[wp.vec3],
    top_idx: wp.array[int],
    top_rest: wp.array[wp.vec3],
    rbar: wp.vec3,
    csum: wp.array[float],
    inv_n: float,
    sums: wp.array[float],  # [sin-weight, cos-weight]
):
    # Best-fit in-plane rotation about +y of the top nodes, taken about the rest centroid
    # (rbar) and the CURRENT centroid (cbar) so the plate is free to translate. Accumulate
    # the cross (sin) and dot (cos) terms in the x-z plane.
    tid = wp.tid()
    p = top_idx[tid]
    cbar_x = csum[0] * inv_n
    cbar_z = csum[2] * inv_n
    ax = top_rest[tid][0] - rbar[0]
    az = top_rest[tid][2] - rbar[2]
    bx = particle_q[p][0] - cbar_x
    bz = particle_q[p][2] - cbar_z
    wp.atomic_add(sums, 0, ax * bz - az * bx)
    wp.atomic_add(sums, 1, ax * bx + az * bz)


@wp.kernel
def compute_theta(sums: wp.array[float], lim: float, relax: float, theta: wp.array[float]):
    # Low-pass the fitted angle toward its new value. The fit<->projection feedback can
    # high-frequency oscillate (the free walls deform, the fit jitters, the projection
    # kicks back); relaxing damps that loop so the hinge angle evolves smoothly.
    fit = wp.clamp(wp.atan2(sums[0], sums[1]), -lim, lim)
    theta[0] = theta[0] + relax * (fit - theta[0])


@wp.kernel
def project_plate(
    top_idx: wp.array[int],
    top_rest: wp.array[wp.vec3],
    rbar: wp.vec3,
    csum: wp.array[float],
    inv_n: float,
    theta: wp.array[float],
    particle_q: wp.array[wp.vec3],
    particle_qd: wp.array[wp.vec3],
):
    # Enforce the rigid plate at the fitted angle + current centroid: rotate (rest-rbar) by
    # theta about +y, then place at the current centroid cbar (which preserves the plate's
    # translation, so projection removes only the non-rigid deformation).
    tid = wp.tid()
    p = top_idx[tid]
    c = wp.cos(theta[0])
    sn = wp.sin(theta[0])
    ax = top_rest[tid][0] - rbar[0]
    ay = top_rest[tid][1] - rbar[1]
    az = top_rest[tid][2] - rbar[2]
    particle_q[p] = wp.vec3(csum[0] * inv_n + ax * c - az * sn,
                            csum[1] * inv_n + ay,
                            csum[2] * inv_n + ax * sn + az * c)
    particle_qd[p] = wp.vec3(0.0, 0.0, 0.0)


@wp.kernel
def zero_array(a: wp.array[float]):
    a[wp.tid()] = 0.0


@wp.kernel
def rotate_hinge_body(theta: wp.array[float], orig_z: float, rbar: wp.vec3,
                      csum: wp.array[float], inv_n: float, idx: int,
                      body_q: wp.array[wp.transform]):
    # Apply the SAME floating-plate transform (rotate the rest origin about rbar by +theta,
    # then translate to the current centroid cbar) to the visual hinge_top body so it tracks
    # the projected top caps and stays coincident with them.
    th = theta[0]
    c = wp.cos(th)
    sn = wp.sin(th)
    ax = 0.0 - rbar[0]
    ay = 0.0 - rbar[1]
    az = orig_z - rbar[2]
    pos = wp.vec3(csum[0] * inv_n + ax * c - az * sn,
                  csum[1] * inv_n + ay,
                  csum[2] * inv_n + ax * sn + az * c)
    body_q[idx] = wp.transform(pos, wp.quat_from_axis_angle(wp.vec3(0.0, 1.0, 0.0), th))


class Example:
    def __init__(self, viewer, args):
        self.viewer = viewer
        self.fps = 60
        self.frame_dt = 1.0 / self.fps
        self.sim_time = 0.0
        self.frame = 0
        self.sim_substeps = int(os.environ.get("SUBSTEPS", "64"))
        self.sim_dt = self.frame_dt / self.sim_substeps
        self.P_MAX = float(os.environ.get("PMAX", "0.5e5"))    # max chamber pressure [Pa] (synthetic mode)
        self.ramp = 90  # frames per actuation phase (synthetic mode)

        # Emergent-bending coupling: the cavity pressure (real, from the experiment log)
        # is the ONLY actuation. The bending angle is integrated from the torque the
        # pressurised chambers exert on the hinge plate via the compliant silicone bond
        # (k_att) -- nothing about the angle is prescribed. gain sets convergence speed
        # (the equilibrium angle is gain-independent); the joint range is clamped.
        # The bending angle emerges from the FEM node displacements (best-fit rotation of
        # the top plate). Only the joint's mechanical range is imposed.
        self.theta_lim = np.radians(float(os.environ.get("THETALIM", "55.0")))  # mechanical range
        self.theta_relax = float(os.environ.get("THETA_RELAX", "0.1"))  # angle low-pass (per substep)
        # Fibre reinforcement (inextensible hoop): resists radial ballooning so pressure
        # goes into axial extension. Off by default -- the explicit penalty destabilises
        # the VBD solve; the free-wall bending already emerges without it.
        self.k_fibre = float(os.environ.get("FIBRE", "0.0"))
        # Compliant base pin: PIN_RELAX<1 under-relaxes the per-substep snap-to-rest so the
        # boundary tet above the bond is never fully sheared (stops the inversion cascade);
        # PIN_VDAMP keeps a fraction of the pinned-node velocity instead of zeroing it.
        self.pin_relax = float(os.environ.get("PIN_RELAX", "1.0"))   # 1.0 = old hard pin
        self.pin_vdamp = float(os.environ.get("PIN_VDAMP", "0.0"))   # 0.0 = old velocity zeroing
        # FREETOP: measure the bend angle but do NOT clamp the top cap (frees chamber length
        # so pressure can extend it instead of buckling the axially over-constrained wall).
        self.freetop = bool(os.environ.get("FREETOP"))
        # Pressure-exclusion band above the bottom bond [m]: faces within this margin are
        # not pressurized (they only shear the pinned-boundary tets). Default 6mm (validated).
        self.press_margin = float(os.environ.get("PRESS_MARGIN_MM", "6.0")) * 0.001

        # Drive the real recorded chamber PRESSURES by default (the experiment's input).
        self.use_data = not os.environ.get("SYNTH")
        self.pscale = float(os.environ.get("PSCALE", "1.0"))
        self.stride = int(os.environ.get("STRIDE", "8"))
        if self.use_data:
            self.d_angle, self.d_p1, self.d_p2 = load_ground_truth(GROUND_TRUTH)
            print(f"driving from ground truth: {len(self.d_angle)} rows | "
                  f"angle 0..{np.degrees(self.d_angle).max():.1f}deg | "
                  f"P1 0..{self.d_p1.max() / 1e5:.2f}bar  P2 0..{self.d_p2.max() / 1e5:.2f}bar")

        mat = json.load(open(MATJSON))
        mu = mat["neo_hookean"]["mu_full"] * 1e6
        # Default k_lambda = real near-incompressible bulk (K - 2/3 mu, K from the fit ~21 MPa
        # -> ~20.8e6). The floating plate lets the chamber extend, so the real stiff volume no
        # longer drives the boundary buckling. NU overrides via lambda=2*mu*nu/(1-2nu); KLAMBDA
        # overrides directly.
        k_bulk = mat["bulk_modulus_MPa"] * 1e6
        nu_env = os.environ.get("NU")
        if nu_env is not None:
            nu = float(nu_env)
            lam_default = 2.0 * mu * nu / (1.0 - 2.0 * nu)
        else:
            lam_default = k_bulk - 2.0 / 3.0 * mu
            nu = (3.0 * k_bulk - 2.0 * mu) / (2.0 * (3.0 * k_bulk + mu))  # for the printout
        k_lambda = float(os.environ.get("KLAMBDA", str(lam_default)))
        # Quasi-static MASS SCALING (gravity off -> static shape is mass-independent):
        # conditions the VBD solve for the tiny tets of this small thin-walled part.
        rho = mat["density_kg_m3"] * float(os.environ.get("DENSCALE", "100.0"))
        print(f"material: k_mu={mu:.3e} k_lambda={k_lambda:.3e} (nu={nu:.3f}) "
              f"tri_ke={float(os.environ.get('TRI', '1.0e5')):.2e} "
              f"denscale={float(os.environ.get('DENSCALE', '100.0')):.0f}")

        tm = newton.TetMesh.create_from_file(MESH)
        tm.custom_attributes = {}
        v = np.asarray(tm.vertices).copy()
        geom = analyze_geometry(v)
        self.geom = geom
        print(f"chamber: {tm.tet_count} tets, {tm.vertex_count} verts | "
              f"L={geom['L']:.1f} r_out={geom['r_out']:.1f} scale={geom['scale']}")
        cav_list = classify_cavities(v, tm.surface_tri_indices, geom)
        cav_local = cav_list[0][0].astype(np.int32)  # the single D-cavity, local indices
        surf_local = np.unique(np.asarray(tm.surface_tri_indices))  # all wall nodes (local)

        builder = newton.ModelBuilder(gravity=0.0)
        builder.particle_max_velocity = 0.3
        s = geom["scale"]
        tri = float(os.environ.get("TRI", "1.0e5"))   # membrane stiffness <= bulk mu (was 8e5)
        kdamp = float(os.environ.get("KDAMP", "2.0e5"))
        # Hoop-spring fibre reinforcement (the inextensible layer).
        self.hoop_ke = float(os.environ.get("HOOP_KE", "2.0e4"))   # hoop spring stiffness [N/m]
        self.hoop_kd = float(os.environ.get("HOOP_KD", "1.0e0"))   # hoop spring damping [N·s/m]
        hoop_band = float(os.environ.get("HOOP_BAND_MM", "4.0")) * s
        hoop_maxlen = float(os.environ.get("HOOP_MAXLEN_MM", "7.0")) * s  # drop cross-cavity springs

        # TWO SEPARATE chambers, one per hinge slot. The +x D-mesh is instanced as-is for
        # the right slot and rotated 180 deg about the length (z) axis for the left slot,
        # so the two bodies sit either side of the central hinge wall with a clear gap
        # between them (no solid silicone bridge across the middle). Mesh is already in
        # the assembly frame (x>=7, z in [0,L]) so no centring offset is applied.
        specs = [("right(+x)", wp.quat_identity()),
                 ("left(-x)", wp.quat_from_axis_angle(wp.vec3(0.0, 0.0, 1.0), float(np.pi)))]
        self.bodies = []
        n_springs = 0
        for name, rot in specs:
            off = builder.particle_count
            builder.add_soft_mesh(
                pos=wp.vec3(0.0, 0.0, 0.0), rot=rot, scale=s, vel=wp.vec3(0.0),
                mesh=tm, density=rho, k_mu=mu, k_lambda=k_lambda, k_damp=kdamp,
                tri_ke=tri, tri_ka=tri, tri_kd=0.05 * tri, validate_mesh=False,
            )
            # Inextensible fibre layer as circumferential hoop springs on this chamber's
            # wall nodes (chamber axis is the cylinder axis through the origin in x-y).
            if self.hoop_ke > 0.0:
                gidx = surf_local + off
                pq = np.asarray(builder.particle_q)[gidx]
                n_springs += build_hoop_springs(builder, gidx, pq, (0.0, 0.0),
                                                hoop_band, self.hoop_ke, self.hoop_kd, hoop_maxlen)
            self.bodies.append({"name": name, "offset": off})
        print(f"  hoop-spring reinforcement: {n_springs} springs (ke={self.hoop_ke:.0f} N/m)")
        # Rigid hinge halves as visual reference (Stage 1). Collision disabled so the
        # tiny-mass chamber particles can't explode against them; gravity is off so the
        # free bodies float in place. Same mm assembly frame -> shape scale s.
        if not os.environ.get("NOHINGE"):
            novis = newton.ModelBuilder.ShapeConfig(has_shape_collision=False, has_particle_collision=False)
            sv = wp.vec3(s, s, s)
            for name, path, col, dz_mm in (
                ("hinge_bottom", HINGE_BOTTOM, wp.vec3(0.30, 0.55, 0.95), 0.0),
                ("hinge_top", HINGE_TOP, wp.vec3(0.95, 0.45, 0.30), HINGE_TOP_DZ_MM),
            ):
                xf = wp.transform(wp.vec3(0.0, 0.0, dz_mm * s), wp.quat_identity())
                b = builder.add_body(xform=xf, label=name)
                builder.add_shape_mesh(b, mesh=load_mesh(path), scale=sv, cfg=novis, color=col, label=name)
            print(f"  loaded rigid hinges (hinge_top raised {HINGE_TOP_DZ_MM}mm; "
                  f"pin axis: +y through x=0, z={PIN_Z_MM}mm)")

        builder.color()
        self.model = builder.finalize()

        rest_q = self.model.particle_q.numpy()
        self.rest_q = rest_q                               # host rest pose for displacement tracking
        self.tet4 = self.model.tet_indices.numpy().reshape(-1, 4)  # [tet_count, 4] connectivity
        self.tet_v0 = self._tet_volumes(rest_q)            # signed rest volumes (all > 0)
        za = rest_q[:, geom["axis"]]
        self.pin_x = 0.0                                   # pin axis x = assembly centre
        self.pin_z = PIN_Z_MM * s                          # z of pin axis (world m)

        # The silicone is bonded to the rigid hinge plates over the chamber ends: nodes
        # BELOW the pin are fixed to hinge_bottom (hard pin = the anchored base); nodes
        # ABOVE the pin are bonded to hinge_top. CRUCIALLY only the chamber END CAPS are
        # bonded to the plates -- the middle wall (z in [z_lo, z_hi]) is left FREE so the
        # cavity pressure can actually inflate/extend it (the silicone is glued to the
        # plates only at its ends, exactly as fabricated). The hinge angle then emerges
        # from how the free walls push the top cap.
        zmax = za.max()
        self.z_lo = float(os.environ.get("ZLO_MM", "8.0")) * s   # bottom bond height
        self.z_hi = zmax - float(os.environ.get("ZHI_MM", "8.0")) * s  # top bond starts here
        bottom = np.where(za < self.z_lo)[0]
        self.pin_idx = wp.array(bottom.astype(np.int32), dtype=int)
        self.pin_rest = wp.array(rest_q[bottom], dtype=wp.vec3)
        self.n_pin = len(bottom)
        self.rest_centroid = rest_q.mean(0)

        top = np.where(za > self.z_hi)[0]
        self.top_idx = wp.array(top.astype(np.int32), dtype=int)
        self.top_rest = wp.array(rest_q[top], dtype=wp.vec3)
        self.n_top = len(top)
        self.n_free = len(rest_q) - self.n_pin - self.n_top
        # Floating-rigid-plate shape matching: rest centroid of the combined top caps, and a
        # device accumulator for the current centroid (recomputed each substep).
        rbar = rest_q[top].mean(0)
        self.top_rbar = wp.vec3(float(rbar[0]), float(rbar[1]), float(rbar[2]))
        self.csum = wp.zeros(3, dtype=float)   # current top-cap centroid accumulator
        self.inv_n_top = 1.0 / float(self.n_top)
        self.theta = wp.zeros(1, dtype=float)  # emergent hinge angle [rad]
        self.sums = wp.zeros(2, dtype=float)   # best-fit rotation accumulators
        self.rest_rad = wp.array(np.hypot(rest_q[:, 0], rest_q[:, 1]).astype(np.float32), dtype=float)
        self.above = wp.array((rest_q[:, 2] > self.pin_z).astype(np.int32), dtype=int)
        # hinge_top is the 2nd rigid body (after hinge_bottom); its rest z origin:
        self.hinge_top_body = 1
        self.hinge_top_z = HINGE_TOP_DZ_MM * s
        self.have_hinge = not os.environ.get("NOHINGE")

        # One pressurized cavity per chamber body (winding baked outward at rest in the
        # world frame), labelled by x-sign so we can drive them antagonistically.
        self.cavities = []
        for b in self.bodies:
            gtris = cav_local + b["offset"]
            centroid = rest_q[np.unique(gtris)].mean(0)
            oriented = bake_outward(gtris, rest_q, centroid)
            side = "A(+x)" if centroid[0] > self.pin_x else "B(-x)"
            v0 = signed_cavity_volume(oriented, rest_q)
            # Pressure-loaded faces EXCLUDE the band just above the bottom bond: those nodes
            # are pinned, so pressure there does no useful work -- it only shears the first
            # free tet layer into inversion. Keep ALL faces for the volume tracker, though.
            triz = rest_q[oriented].mean(axis=1)[:, 2]
            loaded = oriented[triz > self.z_lo + self.press_margin]
            self.cavities.append({
                "q": wp.array(loaded.astype(np.int32), dtype=wp.vec3i),
                "n": len(loaded),
                "p": wp.zeros(1, dtype=float),
                "side": side,
                "tris": oriented.astype(np.int64),  # full cavity, host copy for volume tracker
                "v0": v0,                            # rest cavity volume [m^3] (signed, outward)
            })
        self.cavities.sort(key=lambda c: c["side"])  # A(+x) first, B(-x) second
        print(f"two chambers {[b['name'] for b in self.bodies]} | bottom-cap {self.n_pin} | "
              f"top-cap {self.n_top} (pin z={PIN_Z_MM}mm) | FREE wall {self.n_free} | "
              f"cavities {[(c['side'], c['n']) for c in self.cavities]}")
        # Rest cavity volume: MUST be > 0. A negative value means the pressure normals are
        # baked INWARD, so the follower load sucks the wall in instead of inflating it.
        vol0 = ", ".join(f"{c['side']} {c['v0'] * 1e9:+.0f} mm^3" for c in self.cavities)
        print(f"  rest cavity volume: {vol0}")

        self.solver = newton.solvers.SolverVBD(model=self.model, iterations=int(os.environ.get("ITERS", "24")))
        self.state_0 = self.model.state()
        self.state_1 = self.model.state()
        self.control = self.model.control()
        self.contacts = self.model.contacts()
        self.theta_deg = 0.0
        self.hist_sim = []   # emergent bend angle per frame [deg]
        self.hist_meas = []  # measured ground-truth joint angle per frame [deg]

        size = max(2.0 * geom["r_out"], geom["L"]) * s
        self.viewer.set_model(self.model)
        self.viewer.set_camera(pos=wp.vec3(0.0, -3.0 * size, 0.6 * geom["L"] * s), pitch=-12.0, yaw=90.0)
        self.capture()

    def capture(self):
        if wp.get_device().is_cuda:
            with wp.ScopedCapture() as capture:
                self.simulate()
            self.graph = capture.graph
        else:
            self.graph = None

    def _tet_volumes(self, q):
        """Signed volume of every tet [m^3]. Sign flips (<=0) = an INVERTED/collapsed
        tet -- the 'black hole' that drags neighbours in."""
        a, b, c, d = (q[self.tet4[:, 0]], q[self.tet4[:, 1]],
                      q[self.tet4[:, 2]], q[self.tet4[:, 3]])
        return np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a)) / 6.0

    def diagnostics(self, q, topk=5):
        """Report the worst-displaced nodes and any inverted/collapsed tets so the
        runaway 'black hole' near the pin boundary can be localized."""
        s = self.geom["scale"]
        disp = np.linalg.norm(q - self.rest_q, axis=1)        # per-node displacement [m]
        order = np.argsort(disp)[::-1][:topk]
        zrel = (self.rest_q[:, 2] - self.z_lo) / s            # mm above the bottom bond line
        vol = self._tet_volumes(q)
        ratio = vol / self.tet_v0
        n_inv = int((ratio <= 0.0).sum())                     # inverted tets
        n_crush = int(((ratio > 0.0) & (ratio < 0.1)).sum())  # near-collapsed (<10% rest vol)
        print(f"    DISP[mm] top{topk}: " +
              "  ".join(f"n{n}:{disp[n] * 1e3:.1f}@z{zrel[n]:+.0f}mm" for n in order))
        print(f"    TETS: inverted={n_inv}  crushed(<10%)={n_crush}  "
              f"min vol-ratio={ratio.min():+.2f}")

    def simulate(self):
        for _ in range(self.sim_substeps):
            self.state_0.clear_forces()
            # Cavity pressure (the only actuation) as a real follower load -> the FEM
            # chamber nodes genuinely displace.
            for c in self.cavities:
                wp.launch(apply_cavity_pressure, dim=c["n"],
                          inputs=[self.state_0.particle_q, c["q"], c["p"]],
                          outputs=[self.state_0.particle_f])
            if self.k_fibre > 0.0:
                wp.launch(apply_fibre, dim=self.model.particle_count,
                          inputs=[self.state_0.particle_q, self.above, self.rest_rad,
                                  self.pin_z, self.theta, self.k_fibre],
                          outputs=[self.state_0.particle_f])
            self.model.collide(self.state_0, self.contacts)
            self.solver.step(self.state_0, self.state_1, self.control, self.contacts, self.sim_dt)
            # Fixed base: compliant-pin the bottom nodes to hinge_bottom (under-relaxed so
            # the boundary tet layer cannot be sheared to inversion -- see pin_particles).
            wp.launch(pin_particles, dim=self.n_pin,
                      inputs=[self.pin_idx, self.pin_rest, self.pin_relax, self.pin_vdamp],
                      outputs=[self.state_1.particle_q, self.state_1.particle_qd])
            # Floating rigid plate: fit the best rigid transform (rotation theta + x-z
            # translation) of the displaced top caps, then project them onto it. The
            # translation lets the chambers EXTEND, so the wall no longer buckles; the bend
            # angle theta emerges from the differential extension. (FREETOP skips the
            # projection entirely -- a diagnostic: chambers free, ~0 bend.)
            wp.launch(zero_array, dim=3, inputs=[self.csum])
            wp.launch(accumulate_centroid, dim=self.n_top,
                      inputs=[self.state_1.particle_q, self.top_idx], outputs=[self.csum])
            wp.launch(zero_array, dim=2, inputs=[self.sums])
            wp.launch(accumulate_fit, dim=self.n_top,
                      inputs=[self.state_1.particle_q, self.top_idx, self.top_rest,
                              self.top_rbar, self.csum, self.inv_n_top],
                      outputs=[self.sums])
            wp.launch(compute_theta, dim=1, inputs=[self.sums, self.theta_lim, self.theta_relax],
                      outputs=[self.theta])
            if not self.freetop:
                wp.launch(project_plate, dim=self.n_top,
                          inputs=[self.top_idx, self.top_rest, self.top_rbar, self.csum,
                                  self.inv_n_top, self.theta],
                          outputs=[self.state_1.particle_q, self.state_1.particle_qd])
            if self.have_hinge:
                wp.launch(rotate_hinge_body, dim=1,
                          inputs=[self.theta, self.hinge_top_z, self.top_rbar, self.csum,
                                  self.inv_n_top, self.hinge_top_body],
                          outputs=[self.state_1.body_q])
            self.state_0, self.state_1 = self.state_1, self.state_0

    def schedule(self):
        """Antagonistic demo (normalized actuation per chamber, 0..1):
        ramp A, release, ramp B, release, then co-inflate both."""
        f, R = self.frame, self.ramp
        a = b = 0.0
        if f < R:                      # 1) actuate A  -> bend one way
            a = f / R
        elif f < 2 * R:                # 2) release A
            a = 1 - (f - R) / R
        elif f < 3 * R:                # 3) actuate B  -> bend the other way
            b = (f - 2 * R) / R
        elif f < 4 * R:                # 4) release B
            b = 1 - (f - 3 * R) / R
        else:                          # 5) co-actuate both -> stiffen, ~no net bend
            a = b = min(1.0, (f - 4 * R) / R)
        return a, b

    def step(self):
        if self.use_data:
            i = min(self.frame * self.stride, len(self.d_angle) - 1)
            pA, pB = float(self.d_p1[i]), float(self.d_p2[i])
            self.measured_deg = float(np.degrees(self.d_angle[i]))
            label = f"P1={pA / 1e5:.3f} P2={pB / 1e5:.3f} bar | meas {self.measured_deg:+.1f}"
        else:
            a, b = self.schedule()  # normalized 0..1 per chamber
            pA, pB = a * self.P_MAX, b * self.P_MAX
            self.measured_deg = float("nan")
            label = f"P1={pA / 1e5:.3f} P2={pB / 1e5:.3f} bar"

        # The recorded chamber pressures are the ONLY input. theta is not set here -- it
        # emerges inside simulate() from the pressure-induced torque on the hinge plate.
        self.cavities[0]["p"].fill_(pA * self.pscale)
        self.cavities[1]["p"].fill_(pB * self.pscale)
        if self.graph:
            wp.capture_launch(self.graph)
        else:
            self.simulate()
        self.theta_deg = float(np.degrees(self.theta.numpy()[0]))
        self.hist_sim.append(self.theta_deg)
        self.hist_meas.append(self.measured_deg)
        # Cavity volume tracker: vol/vol0 < 1 means the wall is being pushed INWARD
        # (implosion); > 1 means it is inflating as a real pressurised chamber should.
        q = self.state_0.particle_q.numpy()
        self.vol_ratio = [signed_cavity_volume(c["tris"], q) / c["v0"] for c in self.cavities]
        if self.frame % 10 == 0:
            vr = "  ".join(f"{c['side']} V/V0={r:+.3f}" for c, r in zip(self.cavities, self.vol_ratio))
            print(f"  frame {self.frame:3d}  {label}  theta(sim)={self.theta_deg:+.1f} deg  |  {vr}")
            if os.environ.get("DIAG"):
                self.diagnostics(q)
        self.sim_time += self.frame_dt
        self.frame += 1

    def render(self):
        self.viewer.begin_frame(self.sim_time)
        self.viewer.log_state(self.state_0)
        self.viewer.end_frame()

    def test_final(self):
        q = self.state_0.particle_q.numpy()
        assert np.isfinite(q).all(), "joint diverged (NaN)"
        ext = q.max(0) - q.min(0)
        assert (ext < 0.20).all(), f"joint bbox exploded: {np.round(ext * 1e3, 1)} mm"
        # Ground-truth validation: how well the EMERGENT bend angle tracks the recorded one.
        sim = np.asarray(self.hist_sim)
        meas = np.asarray(self.hist_meas)
        ok = np.isfinite(sim) & np.isfinite(meas)
        if ok.sum() > 2 and np.ptp(meas[ok]) > 1e-6 and np.ptp(sim[ok]) > 1e-9:
            corr = float(np.corrcoef(sim[ok], meas[ok])[0, 1])
            print(f"\n=== bend vs ground truth ===  peak sim {sim[ok].max():+.1f} deg | "
                  f"peak meas {meas[ok].max():+.1f} deg | correlation {corr:+.3f} | "
                  f"sign {'OK' if corr > 0 else 'FLIPPED'}")


if __name__ == "__main__":
    parser = newton.examples.create_parser()
    viewer, args = newton.examples.init(parser)
    newton.examples.run(Example(viewer, args), args)
