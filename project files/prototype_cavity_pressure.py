# SPDX-License-Identifier: Apache-2.0
"""
Pneumatic variable-stiffness joint (Gaozhang 2024) — antagonistic bending.

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
MESH = os.path.join(CADDIR, "chamber_half.msh")
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
    particle_q: wp.array[wp.vec3],
    particle_qd: wp.array[wp.vec3],
):
    # Hard projection pin: particle_inv_mass=0 does NOT hold particles in VBD, so
    # reset pinned nodes to rest pose + zero velocity every substep.
    tid = wp.tid()
    p = pin_idx[tid]
    particle_q[p] = pin_rest[tid]
    particle_qd[p] = wp.vec3(0.0, 0.0, 0.0)


# --- Emergent revolute hinge (pure FEM, angle read from node displacement) ---
# The chambers are FEM; the cavity pressure displaces their nodes. The top of both
# chambers is silicone bonded to a single rigid top plate (hinge_top), which can only
# rotate about the pin axis (+y through x=0, z=pin). Each substep we (1) read the rigid
# rotation that best fits where the pressure pushed the top nodes (a least-squares fit of
# their displacement from rest -- NOT a torque), then (2) project the top nodes onto that
# rigid rotation to enforce the plate. The bending angle therefore EMERGES from the FEM
# node displacements; the silicone spanning the pin supplies the real restoring, so the
# joint settles at the physical equilibrium with nothing about the angle prescribed.
@wp.func
def _plate_target(rest: wp.vec3, pin_x: float, pin_z: float, c: float, sn: float) -> wp.vec3:
    rx = rest[0] - pin_x
    rz = rest[2] - pin_z
    return wp.vec3(pin_x + rx * c - rz * sn, rest[1], pin_z + rx * sn + rz * c)


@wp.kernel
def accumulate_rotation(
    particle_q: wp.array[wp.vec3],
    top_idx: wp.array[int],
    top_rest: wp.array[wp.vec3],
    pin_x: float,
    pin_z: float,
    sums: wp.array[float],  # [sin-weight, cos-weight]
):
    # Best-fit rotation about +y of the displaced top nodes relative to their rest pose:
    # accumulate the cross (sin) and dot (cos) terms in the x-z plane about the pin.
    tid = wp.tid()
    p = top_idx[tid]
    rx = top_rest[tid][0] - pin_x
    rz = top_rest[tid][2] - pin_z
    cx = particle_q[p][0] - pin_x
    cz = particle_q[p][2] - pin_z
    wp.atomic_add(sums, 0, rx * cz - rz * cx)
    wp.atomic_add(sums, 1, rx * cx + rz * cz)


@wp.kernel
def compute_theta(sums: wp.array[float], lim: float, theta: wp.array[float]):
    theta[0] = wp.clamp(wp.atan2(sums[0], sums[1]), -lim, lim)


@wp.kernel
def project_top(
    top_idx: wp.array[int],
    top_rest: wp.array[wp.vec3],
    pin_x: float,
    pin_z: float,
    theta: wp.array[float],
    particle_q: wp.array[wp.vec3],
    particle_qd: wp.array[wp.vec3],
):
    # Enforce the rigid top plate at the emergent (fitted) angle.
    tid = wp.tid()
    p = top_idx[tid]
    target = _plate_target(top_rest[tid], pin_x, pin_z, wp.cos(theta[0]), wp.sin(theta[0]))
    particle_q[p] = target
    particle_qd[p] = wp.vec3(0.0, 0.0, 0.0)


@wp.kernel
def zero_array(a: wp.array[float]):
    a[wp.tid()] = 0.0


@wp.kernel
def rotate_hinge_body(theta: wp.array[float], orig_z: float, pin_z: float, idx: int,
                      body_q: wp.array[wp.transform]):
    # Rotate hinge_top about the SAME pin axis the chamber top is projected about, so the
    # plate and the chamber top cap share one rotation centre and stay coincident.
    # project_top maps rest -> R_y(-theta)*(rest-pin)+pin, so the plate must use the
    # SAME rotation R_y(-theta) (a +theta quat would spin it the opposite way).
    th = theta[0]
    c = wp.cos(th)
    sn = wp.sin(th)
    rz = orig_z - pin_z
    body_q[idx] = wp.transform(wp.vec3(-rz * sn, 0.0, pin_z + rz * c),
                               wp.quat_from_axis_angle(wp.vec3(0.0, 1.0, 0.0), -th))


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
        # Fibre reinforcement (inextensible hoop): resists radial ballooning so pressure
        # goes into axial extension. Off by default -- the explicit penalty destabilises
        # the VBD solve; the free-wall bending already emerges without it.
        self.k_fibre = float(os.environ.get("FIBRE", "0.0"))

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
        nu = 0.45
        k_lambda = 2.0 * mu * nu / (1.0 - 2.0 * nu)
        # Quasi-static MASS SCALING (gravity off -> static shape is mass-independent):
        # conditions the VBD solve for the tiny tets of this small thin-walled part.
        rho = mat["density_kg_m3"] * float(os.environ.get("DENSCALE", "50.0"))

        tm = newton.TetMesh.create_from_file(MESH)
        tm.custom_attributes = {}
        v = np.asarray(tm.vertices).copy()
        geom = analyze_geometry(v)
        self.geom = geom
        print(f"chamber: {tm.tet_count} tets, {tm.vertex_count} verts | "
              f"L={geom['L']:.1f} r_out={geom['r_out']:.1f} scale={geom['scale']}")
        cav_list = classify_cavities(v, tm.surface_tri_indices, geom)
        cav_local = cav_list[0][0].astype(np.int32)  # the single D-cavity, local indices

        builder = newton.ModelBuilder(gravity=0.0)
        builder.particle_max_velocity = 0.3
        s = geom["scale"]
        tri = float(os.environ.get("TRI", "8.0e5"))
        kdamp = float(os.environ.get("KDAMP", "1.0e5"))

        # TWO SEPARATE chambers, one per hinge slot. The +x D-mesh is instanced as-is for
        # the right slot and rotated 180 deg about the length (z) axis for the left slot,
        # so the two bodies sit either side of the central hinge wall with a clear gap
        # between them (no solid silicone bridge across the middle). Mesh is already in
        # the assembly frame (x>=7, z in [0,L]) so no centring offset is applied.
        specs = [("right(+x)", wp.quat_identity()),
                 ("left(-x)", wp.quat_from_axis_angle(wp.vec3(0.0, 0.0, 1.0), float(np.pi)))]
        self.bodies = []
        for name, rot in specs:
            off = builder.particle_count
            builder.add_soft_mesh(
                pos=wp.vec3(0.0, 0.0, 0.0), rot=rot, scale=s, vel=wp.vec3(0.0),
                mesh=tm, density=rho, k_mu=mu, k_lambda=k_lambda, k_damp=kdamp,
                tri_ke=tri, tri_ka=tri, tri_kd=0.05 * tri, validate_mesh=False,
            )
            self.bodies.append({"name": name, "offset": off})
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
            self.cavities.append({
                "q": wp.array(oriented.astype(np.int32), dtype=wp.vec3i),
                "n": len(oriented),
                "p": wp.zeros(1, dtype=float),
                "side": side,
            })
        self.cavities.sort(key=lambda c: c["side"])  # A(+x) first, B(-x) second
        print(f"two chambers {[b['name'] for b in self.bodies]} | bottom-cap {self.n_pin} | "
              f"top-cap {self.n_top} (pin z={PIN_Z_MM}mm) | FREE wall {self.n_free} | "
              f"cavities {[(c['side'], c['n']) for c in self.cavities]}")

        self.solver = newton.solvers.SolverVBD(model=self.model, iterations=int(os.environ.get("ITERS", "20")))
        self.state_0 = self.model.state()
        self.state_1 = self.model.state()
        self.control = self.model.control()
        self.contacts = self.model.contacts()
        self.theta_deg = 0.0

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
            # Fixed base: hard-pin the bottom nodes to hinge_bottom.
            wp.launch(pin_particles, dim=self.n_pin,
                      inputs=[self.pin_idx, self.pin_rest],
                      outputs=[self.state_1.particle_q, self.state_1.particle_qd])
            # Emergent angle: fit the rigid rotation of the displaced top nodes, then
            # project them onto that rotation to enforce the shared rigid top plate.
            wp.launch(zero_array, dim=2, inputs=[self.sums])
            wp.launch(accumulate_rotation, dim=self.n_top,
                      inputs=[self.state_1.particle_q, self.top_idx, self.top_rest, self.pin_x, self.pin_z],
                      outputs=[self.sums])
            wp.launch(compute_theta, dim=1, inputs=[self.sums, self.theta_lim], outputs=[self.theta])
            wp.launch(project_top, dim=self.n_top,
                      inputs=[self.top_idx, self.top_rest, self.pin_x, self.pin_z, self.theta],
                      outputs=[self.state_1.particle_q, self.state_1.particle_qd])
            if self.have_hinge:
                wp.launch(rotate_hinge_body, dim=1,
                          inputs=[self.theta, self.hinge_top_z, self.pin_z, self.hinge_top_body],
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
        if self.frame % 30 == 0:
            print(f"  frame {self.frame:3d}  {label}  theta(sim)={self.theta_deg:+.1f} deg")
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


if __name__ == "__main__":
    parser = newton.examples.create_parser()
    viewer, args = newton.examples.init(parser)
    newton.examples.run(Example(viewer, args), args)
