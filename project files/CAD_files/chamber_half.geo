// ============================================================
// SINGLE semicircular (D-shaped) actuation chamber for one slot
// of the Gaozhang-2024 joint. This is instanced TWICE in the sim
// (right +x as-is, left -x rotated 180 deg about the length axis)
// to form the two SEPARATE antagonistic chambers either side of the
// central rigid hinge wall -- they do NOT join across the middle.
//
//   outer arc radius R   = 33 mm   (2 mm inside the 35 mm flange)
//   inner cavity arc Rin = 29 mm   (4 mm outer/arc wall)
//   slot gap             = x > 7 mm (clears central hinge wall |x|<6.25)
//   chord wall           = 4 mm    (cavity flat face at x = 11 mm)
//   length L             = 69 mm   (z 0 -> 69, between the two flanges)
//   end caps             = 4 mm    (cavity spans z in [4, 65])
//   => pin axis (z=34.6 mm) sits at chamber mid-height
// One closed semicircular cavity (D cross-section). Inflating it
// extends that slot axially -> the joint bends about the pin.
//   gmsh chamber_half.geo -3 -format msh2 -o chamber_half.msh
// ============================================================
SetFactory("OpenCASCADE");

R    = 33.0;
Rin  = 29.0;
gap  = 7.0;     // inner flat face of the body (clears the hinge wall)
wall = 4.0;     // chord-wall thickness  -> cavity flat face at gap+wall
L    = 69.0;
cap  = 4.0;
lc   = 3.5;
eps  = 1e-3;

Mesh.CharacteristicLengthMin = lc * 0.6;
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
