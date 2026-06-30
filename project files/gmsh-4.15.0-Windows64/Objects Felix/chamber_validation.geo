SetFactory("OpenCASCADE");

// ---------------- Chamber Parameters ----------------
// Matching our soft joint chamber: 38mm diameter, 49mm height, 2.5mm wall
R  = 19;     // outer radius (mm) — half of 38mm diameter
L  = 49;     // height (mm)
t  = 2.5;    // wall thickness (mm) → inner radius = 16.5mm
lc = 3.0;    // mesh size (mm) — coarser than Felix to match our ~400 element target
eps = 1e-3;

Mesh.CharacteristicLengthMin = lc * 0.5;
Mesh.CharacteristicLengthMax = lc * 1.5;

// ---------------- 1) Outer half-cylinder ----------------
Cylinder(1) = {0, 0, 0,  0, 0, L,  R};
// Keep the half where x >= 0 (right half, our chamber orientation)
Box(2) = {-eps, -R-eps, -eps,  R+2*eps, 2*(R+eps), L+2*eps};
outerHalf[] = BooleanIntersection{ Volume{1}; Delete; }{ Volume{2}; Delete; };

// ---------------- 2) Inner cavity (same geometry, inset by t) ----------------
Rin = R - t;  // 16.5 mm
Lin = L;      // full height (cavity open at top/bottom)

Cylinder(3) = {0, 0, 0,  0, 0, Lin,  Rin};
Box(4) = {-eps, -Rin-eps, -eps,  Rin+2*eps, 2*(Rin+eps), Lin+2*eps};
innerHalf[] = BooleanIntersection{ Volume{3}; Delete; }{ Volume{4}; Delete; };

// ---------------- 3) Subtract cavity from outer solid ----------------
solid[] = BooleanDifference{ Volume{outerHalf[]}; Delete; }{ Volume{innerHalf[]}; Delete; };

// ---------------- 4) Mesh ----------------
Physical Volume("Chamber") = {solid[]};

// Tetrahedral mesh settings
Mesh.Algorithm3D = 1;      // 1 = Delaunay (good quality for thin walls)
Mesh.Recombine3DAll = 0;
Mesh.RecombineAll = 0;
Mesh.ElementOrder = 1;     // linear tets (Tet4)
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Mesh 3;
Mesh.SaveAll = 1;
