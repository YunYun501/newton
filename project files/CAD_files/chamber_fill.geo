// ============================================================
// Joint chamber sized to FILL the user's hinge assembly.
//   outer radius R   = 33 mm   (just inside the 35 mm flange)
//   cavity radius Rin= 29 mm    (4 mm outer wall)
//   central web      = +/- 7 mm (clears the hinge knuckles |x|<7)
//   length L         = 69 mm    (z 0 -> 69, bottom flange top to top flange)
//   end caps         = 4 mm
//   => pin axis (z=34.6 mm) sits at the chamber mid-height
// Two semicircular cavities (+x / -x) separated by the central web.
//   gmsh chamber_fill.geo -3 -format msh2 -o chamber_fill.msh
// ============================================================
SetFactory("OpenCASCADE");

R   = 33.0;
Rin = 29.0;
web = 7.0;
L   = 69.0;
cap = 4.0;
lc  = 3.5;
eps = 1e-3;

Mesh.CharacteristicLengthMin = lc * 0.6;
Mesh.CharacteristicLengthMax = lc;

Cylinder(1) = {0, 0, 0,  0, 0, L,  R};

Cylinder(2) = {0, 0, cap,  0, 0, L-2*cap,  Rin};
Box(3) = {web, -Rin-eps, cap-eps,  R, 2*(Rin+eps), L-2*cap+2*eps};
BooleanIntersection(4) = { Volume{2}; Delete; }{ Volume{3}; Delete; };

Cylinder(5) = {0, 0, cap,  0, 0, L-2*cap,  Rin};
Box(6) = {-web-R, -Rin-eps, cap-eps,  R, 2*(Rin+eps), L-2*cap+2*eps};
BooleanIntersection(7) = { Volume{5}; Delete; }{ Volume{6}; Delete; };

BooleanDifference(8) = { Volume{1}; Delete; }{ Volume{4}; Volume{7}; Delete; };
Physical Volume("Joint") = {8};

Mesh.Algorithm3D = 1;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.ElementOrder = 1;
Mesh 3;
Mesh.SaveAll = 1;
