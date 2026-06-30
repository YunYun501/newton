// ============================================================
// Gmsh .geo — Hollow Half-Cylinder Chamber
//
// Reconstructed parametrically from Chamber_STL.STL geometry:
//   Circular arch cross-section (~157 deg arc, not full semicircle)
//   Outer radius: 35    Inner radius: 31    Wall: 4
//   Circle center: (34.2929, 35.0)
//   Outer chord (flat top): Y = 28      Inner chord: Y = 24
//   Total Z length: 75    Solid end caps: 4 thick
//
// Usage:
//   Command line:  gmsh Chamber_mesh.geo -3
//   GUI:           Open in Gmsh -> Mesh -> 3D
// ============================================================
SetFactory("OpenCASCADE");

// --- Geometry parameters (edit these to resize) ---
R_out = 35.0;           // Outer radius
R_in  = 31.0;           // Inner radius  (wall = R_out - R_in = 4)
cx    = 34.2929;        // Circle center X
cy    = 35.0;           // Circle center Y
L     = 75.0;           // Total length along Z
cap   = 4.0;            // End-cap thickness

Y_chord_out = 28.0;     // Outer arch chord (flat top)
Y_chord_in  = 24.0;     // Inner arch chord

// --- Mesh sizing ---
lc = 10;   // Target element size (~2 elements through wall)

// =============================================================
// Step 1: Outer arch solid (Z = 0 to 75)
// =============================================================
Cylinder(1) = {cx, cy, 0,   0, 0, L,   R_out};

Box(2) = {cx - R_out - 10, Y_chord_out, -10,
           2*R_out + 20,   R_out + 20,   L + 20};

BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

// =============================================================
// Step 2: Inner cavity (Z = cap to L-cap)
// =============================================================
Cylinder(4) = {cx, cy, cap,   0, 0, L - 2*cap,   R_in};

Box(5) = {cx - R_in - 10, Y_chord_in, cap - 10,
           2*R_in + 20,   R_in + 20,   L - 2*cap + 20};

BooleanDifference(6) = { Volume{4}; Delete; }{ Volume{5}; Delete; };

// =============================================================
// Step 3: Subtract cavity from solid -> hollow chamber wall
// =============================================================
BooleanDifference(7) = { Volume{3}; Delete; }{ Volume{6}; Delete; };

// =============================================================
// Physical groups (for SOFA / FEM export)
// =============================================================
Physical Volume("ChamberWall") = {7};
Physical Surface("ChamberSurface") = Surface{:};

// =============================================================
// Mesh settings
// =============================================================
Mesh.CharacteristicLengthMin = lc * 0.5;
Mesh.CharacteristicLengthMax = lc;
Mesh.Algorithm   = 6;    // Frontal-Delaunay (2D)
Mesh.Algorithm3D = 1;    // Delaunay (3D)
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.ElementOrder = 1;   // 1 = linear tets, 2 = quadratic

// Uncomment if SOFA needs older format:
// Mesh.MshFileVersion = 2.2;
