// ============================================================
// Paper-accurate variable-stiffness joint (Gaozhang 2024)
// ONE silicone body with TWO semicircular pneumatic cavities
// separated by a central solid web (the hinge region).
//   outer radius R   = 19 mm  (38 mm diameter)
//   cavity radius Rin= 16.5 mm (2.5 mm outer wall)
//   central web      = +/- 2 mm (4 mm solid strip at x=0 -> bending axis)
//   length L         = 49 mm
//   end caps         = 3 mm    -> cavities span z in [3, 46]
// Cavity 1 lives in the +x half, cavity 2 in the -x half. Inflating one
// extends that side axially -> the body bends about the central web.
//   gmsh chamber_joint.geo -3 -format msh2 -o chamber_joint.msh
// ============================================================
SetFactory("OpenCASCADE");

R   = 19.0;
Rin = 16.5;
web = 2.0;        // half-width of the central solid web
L   = 49.0;
cap = 3.0;
lc  = 2.5;
eps = 1e-3;

Mesh.CharacteristicLengthMin = lc * 0.6;
Mesh.CharacteristicLengthMax = lc;

// ---- full outer body ----
Cylinder(1) = {0, 0, 0,  0, 0, L,  R};

// ---- cavity 1 (+x half): r < Rin, x > web, z in [cap, L-cap] ----
Cylinder(2) = {0, 0, cap,  0, 0, L-2*cap,  Rin};
Box(3) = {web, -Rin-eps, cap-eps,  R, 2*(Rin+eps), L-2*cap+2*eps};
BooleanIntersection(4) = { Volume{2}; Delete; }{ Volume{3}; Delete; };

// ---- cavity 2 (-x half): r < Rin, x < -web, z in [cap, L-cap] ----
Cylinder(5) = {0, 0, cap,  0, 0, L-2*cap,  Rin};
Box(6) = {-web-R, -Rin-eps, cap-eps,  R, 2*(Rin+eps), L-2*cap+2*eps};
BooleanIntersection(7) = { Volume{5}; Delete; }{ Volume{6}; Delete; };

// ---- joint body = full cylinder minus both cavities ----
BooleanDifference(8) = { Volume{1}; Delete; }{ Volume{4}; Volume{7}; Delete; };
Physical Volume("Joint") = {8};

Mesh.Algorithm3D = 1;     // Delaunay
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.ElementOrder = 1;
Mesh 3;
Mesh.SaveAll = 1;
