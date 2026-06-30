// ============================================================
// Gmsh .geo — Hollow Half-Cylinder Chamber (RIGHT / Mirrored)
//
// Mirror of Chamber_mesh.geo across the flat chord plane Y = 28
//   Original arch: Y = 0  → 28  (downward)
//   Mirrored arch: Y = 28 → 56  (upward)
//
//   Circle center mirrored: (34.2929, 21.0)  [was (34.2929, 35.0)]
//   Outer radius: 35    Inner radius: 31    Wall: 4
//   Outer chord: Y = 28    Inner chord: Y = 32
//   Total Z length: 75    Solid end caps: 4 thick
//
// Usage:
//   gmsh Chamber_Right.geo -3
// ============================================================
SetFactory("OpenCASCADE");

// --- Geometry parameters ---
R_out = 35.0;
R_in  = 31.0;
cx    = 34.2929;
cy    = 21.0;           // Mirrored center Y  (original 35, mirror across Y=28: 2*28-35 = 21)
L     = 75.0;
cap   = 4.0;

Y_chord_out = 28.0;     // Outer chord (shared mirror plane)
Y_chord_in  = 32.0;     // Inner chord (mirrored from Y=24: 2*28-24 = 32)

// --- Mesh sizing ---
lc = 10;

// =============================================================
// Step 1: Outer arch solid (Z = 0 to 75)
// =============================================================
Cylinder(1) = {cx, cy, 0,   0, 0, L,   R_out};

// Cut box removes everything BELOW the chord (Y < 28)
Box(2) = {cx - R_out - 10, cy - R_out - 10, -10,
           2*R_out + 20,   Y_chord_out - (cy - R_out - 10),   L + 20};

BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };
// Result: upward arch, Y in [28, 56]

// =============================================================
// Step 2: Inner cavity (Z = cap to L-cap)
// =============================================================
Cylinder(4) = {cx, cy, cap,   0, 0, L - 2*cap,   R_in};

// Cut box removes everything BELOW inner chord (Y < 32)
Box(5) = {cx - R_in - 10, cy - R_in - 10, cap - 10,
           2*R_in + 20,   Y_chord_in - (cy - R_in - 10),   L - 2*cap + 20};

BooleanDifference(6) = { Volume{4}; Delete; }{ Volume{5}; Delete; };
// Result: inner upward arch, Y in [32, 52]

// =============================================================
// Step 3: Subtract cavity -> hollow chamber wall
// =============================================================
BooleanDifference(7) = { Volume{3}; Delete; }{ Volume{6}; Delete; };

// =============================================================
// Physical groups
// =============================================================
Physical Volume("ChamberWall") = {7};
Physical Surface("ChamberSurface") = Surface{:};

// =============================================================
// Mesh settings
// =============================================================
Mesh.CharacteristicLengthMin = lc * 0.5;
Mesh.CharacteristicLengthMax = lc;
Mesh.Algorithm   = 6;
Mesh.Algorithm3D = 1;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.ElementOrder = 1;
