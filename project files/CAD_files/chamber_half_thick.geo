// ============================================================
// Thick-wall D-shaped chamber for the hoop-spring reinforced joint.
// Thicker walls + finer mesh give the wall real BENDING stiffness and
// >=2-3 tet layers through the thickness, so it resists buckling at the
// full recorded pressure (~0.74 bar) when paired with the hoop springs.
//
//   outer arc radius R   = 33 mm   (2 mm inside the 35 mm flange)
//   inner cavity arc Rin = 27 mm   (6 mm outer/arc wall, was 4 mm)
//   slot gap             = x > 7 mm (clears central hinge wall |x|<6.25)
//   chord wall           = 6 mm    (cavity flat face at x = 13 mm)
//   length L             = 69 mm   (z 0 -> 69, between the two flanges)
//   end caps             = 6 mm    (cavity spans z in [6, 63])
//   lc                   = 2.5 mm  (~2-3 tets across the 6 mm wall)
//   gmsh chamber_half_thick.geo -3 -format msh2 -o chamber_half_thick.msh
// ============================================================
SetFactory("OpenCASCADE");

R    = 33.0;
Rin  = 27.0;
gap  = 7.0;
wall = 6.0;
L    = 69.0;
cap  = 6.0;
lc   = 2.5;
eps  = 1e-3;

Mesh.CharacteristicLengthMin = lc * 0.5;
Mesh.CharacteristicLengthMax = lc;

// ---- outer D-shaped solid: cylinder(R) intersect {x > gap} ----
Cylinder(1) = {0, 0, 0,  0, 0, L,  R};
Box(2) = {gap, -R-eps, -eps,  R, 2*(R+eps), L+2*eps};
BooleanIntersection(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

// ---- closed cavity: cylinder(Rin) intersect {x > gap+wall} intersect {z in [cap,L-cap]} ----
Cylinder(4) = {0, 0, cap,  0, 0, L-2*cap,  Rin};
Box(5) = {gap+wall, -Rin-eps, cap-eps,  Rin, 2*(Rin+eps), L-2*cap+2*eps};
BooleanIntersection(6) = { Volume{4}; Delete; }{ Volume{5}; Delete; };

// ---- chamber = outer solid minus cavity ----
BooleanDifference(7) = { Volume{3}; Delete; }{ Volume{6}; Delete; };
Physical Volume("Chamber") = {7};

Mesh.Algorithm3D = 1;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.ElementOrder = 1;
Mesh 3;
Mesh.SaveAll = 1;
