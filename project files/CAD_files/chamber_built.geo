// ============================================================
// Paper-accurate semicircular pneumatic chamber (Gaozhang 2024)
// Closed hollow half-cylinder with solid end caps.
//   outer radius R = 19 mm  (38 mm diameter)
//   wall t       = 2.5 mm   -> inner radius Rin = 16.5 mm
//   length L     = 49 mm
//   end caps     = 3 mm     -> cavity spans z in [3, 46]
// Flat (chord) side at x = 0; cavity inset by the wall on every side
// so the air void is fully enclosed (a closed component for the
// connected-component cavity classifier).
//   gmsh chamber_built.geo -3 -o chamber_built.msh
// ============================================================
SetFactory("OpenCASCADE");

R   = 19.0;
t   = 2.5;
Rin = R - t;      // 16.5
L   = 49.0;
cap = 3.0;
lc  = 2.0;        // element size (mm): ~2 layers across the 2.5 mm wall
eps = 1e-3;

Mesh.CharacteristicLengthMin = lc * 0.6;
Mesh.CharacteristicLengthMax = lc;

// ---- outer half-cylinder (keep x >= 0) ----
Cylinder(1) = {0, 0, 0,  0, 0, L,  R};
Box(2) = {-eps, -R-eps, -eps,  R+2*eps, 2*(R+eps), L+2*eps};
BooleanIntersection(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

// ---- inner cavity half-cylinder (keep x >= t), z in [cap, L-cap] ----
Cylinder(4) = {0, 0, cap,  0, 0, L-2*cap,  Rin};
Box(5) = {t, -Rin-eps, cap-eps,  R+2*eps, 2*(Rin+eps), L-2*cap+2*eps};
BooleanIntersection(6) = { Volume{4}; Delete; }{ Volume{5}; Delete; };

// ---- hollow chamber = outer solid minus cavity ----
BooleanDifference(7) = { Volume{3}; Delete; }{ Volume{6}; Delete; };
Physical Volume("Chamber") = {7};

Mesh.Algorithm3D = 1;     // Delaunay (good for thin walls)
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.ElementOrder = 1;    // linear tets
Mesh 3;
Mesh.SaveAll = 1;
