# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

@AGENTS.md

The development conventions, lint/format/test commands, and PR rules live in `AGENTS.md`
(imported above). This file adds the high-level architecture that those rules assume.

## What Newton is

Newton is a GPU-accelerated physics engine for robotics, built on **NVIDIA Warp**.
Nearly all compute happens in Warp kernels (`@wp.kernel` / `@wp.func`) operating on
`wp.array` data that lives on-device (CPU or CUDA). Python code is mostly orchestration:
it builds data, launches kernels, and swaps buffers. When reading core code, expect the
"real" logic to be in kernels, with the surrounding Python class being a thin owner of
arrays. MuJoCo Warp (`mujoco_warp`) is integrated as the primary rigid-body backend.

## Public API vs. internal `_src`

The single most important structural rule (enforced by ruff `TID253` and `test_api.py`):

- **`newton/_src/`** holds all implementation. Nothing user-facing should be imported from here.
- **`newton/*.py`** are the public façade modules. Each re-exports a curated set of symbols
  from `_src` and defines `__all__`. The mapping is direct:
  - `newton/__init__.py` → core types + the central sim objects (`Model`, `ModelBuilder`,
    `State`, `Control`, `Contacts`, joint/flag enums, `eval_fk`/`eval_ik`/...).
  - `newton/solvers.py` → `_src/solvers` (the `SolverBase` subclasses).
  - `newton/geometry.py`, `actuators.py`, `ik.py`, `math.py`, `selection.py`,
    `sensors.py`, `usd.py`, `utils.py`, `viewer.py` → their matching `_src/<name>` package.

When adding a public symbol: implement in `_src`, re-export in the matching public module,
add it to `__all__`, then run `docs/generate_api.py`. Examples, docs, and tests import only
the public modules. Heavy/optional third-party imports (mujoco, pxr, trimesh, torch, pyglet,
...) are banned at module level in `_src` and `tests` (ruff `TID253`) — import them lazily
inside functions.

## The simulation model: build → finalize → step

The core data flow is the same across every example and solver:

1. **`ModelBuilder`** (`_src/sim/builder.py`) accumulates a scene on the host: bodies/links,
   shapes (`add_shape_box`, `add_shape_sphere`, ...), joints (`add_joint_revolute`, ...),
   particles, articulations, and ground plane. Importers (`_src/utils/import_{urdf,mjcf,usd}.py`)
   populate a builder from asset files.
2. **`builder.finalize()`** produces an immutable **`Model`** (`_src/sim/model.py`) — the
   on-device arrays describing the scene (masses, inertias, joint params, shape geometry).
   Pass `requires_grad=True` here to enable differentiable rollouts.
3. **`Model.state()`** / **`.control()`** / **`.contacts()`** mint the mutable per-step buffers:
   - **`State`** (`_src/sim/state.py`): positions/velocities/forces — `joint_q`, `joint_qd`,
     `body_q`, `body_qd`, `particle_q`, etc. This is what advances in time.
   - **`Control`** (`_src/sim/control.py`): actuation inputs (joint targets, feedforward forces).
   - **`Contacts`** (`_src/sim/contacts.py`): contact constraints produced by collision.
4. The **step loop** (see `newton/examples/basic/example_basic_pendulum.py` for the canonical
   shape):
   ```python
   state_0.clear_forces()
   model.collide(state_0, contacts)                       # broad+narrow phase → Contacts
   solver.step(state_0, state_1, control, contacts, dt)   # integrate one substep
   state_0, state_1 = state_1, state_0                    # double-buffer swap
   ```
   On CUDA the substep loop is captured once into a `wp.ScopedCapture` graph and replayed
   per frame for performance.

Coordinate conventions matter: articulation solvers (`SolverFeatherstone`, `SolverMuJoCo`)
use **generalized/reduced** coordinates; the others (`SolverSemiImplicit`, `SolverXPBD`,
`SolverKamino`) use **maximal** coordinates and enforce joints as pairwise constraints.
`eval_fk` populates `body_q`/`body_qd` from `joint_q`/`joint_qd` (needed by non-MuJoCo
solvers after construction).

## Solvers

All solvers subclass **`SolverBase`** (`_src/solvers/solver.py`) and implement `step()`.
Each lives in its own `_src/solvers/<name>/solver_<name>.py`:

- `SolverMuJoCo` — primary rigid-body / articulation backend via `mujoco_warp`.
- `SolverFeatherstone` — reduced-coordinate articulation dynamics.
- `SolverSemiImplicit`, `SolverXPBD` — maximal-coordinate, support particles/cloth/soft bodies.
- `SolverVBD`, `SolverStyle3D` — implicit, cloth/soft-body focused.
- `SolverImplicitMPM` — material point method for granular/continuum media.
- `SolverKamino` — experimental maximal-coordinate solver (has its own nested `_src` tree).

`newton/solvers.py`'s module docstring holds the authoritative feature-support matrix
(which solver supports which joint type / property / integration scheme). Update it when
changing solver capabilities. `SolverNotifyFlags` tells a solver which `Model` arrays changed
so it can refresh cached state without a full rebuild.

## Other `_src` subsystems

- **`geometry/`** — shape types (`Mesh`, `SDF`, `Heightfield`, ...), the collision pipeline
  (broad phase `broad_phase_{nxn,sap}.py`, narrow phase, GJK/MPR, contact reduction,
  hydroelastic/SDF contacts), and inertia computation.
- **`actuators/`** — actuator library (`ActuatorPD`, controllers, clamping models) layered on
  top of `Control`.
- **`sim/ik/`** + `newton/ik.py` — inverse kinematics objectives and LM/L-BFGS optimizers.
- **`sensors/`** — contact, IMU, frame-transform, and tiled-camera sensors (the camera path
  has its own Warp ray-tracer under `sensors/warp_raytrace/`).
- **`utils/`** — asset importers, mesh processing, asset download, rendering helpers.
- **`viewer/`** + `newton/viewer.py` — viewers selected by the `--viewer` flag
  (`gl`, `usd`, `rtx`, `rerun`, `viser`, `null`).

## Examples

Examples live under `newton/examples/<category>/example_<name>.py` and are run via
`python -m newton.examples <name>` (see `README.md` for the full catalog and `--list`).
Each follows the **`Example` class** contract: `__init__(viewer, args)` builds the model,
`step()` advances one frame, `render()` logs to the viewer, and `test_final()` (plus optional
`test_post_step()`) validate simulation state — these are what `tests/test_examples.py` runs in
CI. Use `newton.examples.test_body_state(...)` for state assertions. Register new examples in
`README.md` with a command and a 320×320 jpg.

## Tests

Tests use **`unittest`** (never pytest), run via `python -m newton.tests` (a parallel runner).
`tests/test_api.py` guards the public surface; `tests/test_generate_api.py` checks the docs
generator is in sync. See `AGENTS.md` for the exact `uv run` invocations, including the
single-test `-k` form and the GPU/torch extras.
