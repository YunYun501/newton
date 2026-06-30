SetFactory("OpenCASCADE");

// ---------------- Parameters ----------------
R  = 10;     // outer radius
L  = 30;     // length along z
t  = 2;      // wall thickness (offset distance)
lc = 1.0;    // mesh size
eps = 1e-3;

// Constraints for a non-empty cavity:
// 1) R > 2*t   (otherwise inner half cross-section vanishes)
// 2) L > 2*t   (otherwise no cavity in z)

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// ---------------- 1) Outer half-cylinder (the outside shape) ----------------
Cylinder(1) = {0, 0, 0,  0, 0, L,  R};

// keep y >= 0 half
Box(2) = {-R-eps, 0, -eps,  2*(R+eps), (R+eps), L+2*eps};
outerHalf[] = BooleanIntersection{ Volume{1}; Delete; }{ Volume{2}; Delete; };

// ---------------- 2) Inner "offset" cavity (same geometry, inset by t) ----------------
Rin = R - t;
Lin = L - 2*t;

// inner cylinder is shorter in z (top/bottom offset by t)
Cylinder(3) = {0, 0, t,  0, 0, Lin,  Rin};

// inner half is also offset from the cut plane by t: keep y >= t
Box(4) = {-Rin-eps, t, t-eps,  2*(Rin+eps), (Rin - t + 2*eps), Lin+2*eps};
innerHalf[] = BooleanIntersection{ Volume{3}; Delete; }{ Volume{4}; Delete; };

// ---------------- 3) Subtract cavity from outer solid ----------------
solid[] = BooleanDifference{ Volume{outerHalf[]}; Delete; }{ Volume{innerHalf[]}; Delete; };

// ---------------- 4) Mesh ----------------
Physical Volume("HalfShellSolid") = {solid[]};

// Force tetrahedral volume meshing
Mesh.Algorithm3D = 10;     // 10 = HXT (tetra), fast/parallel :contentReference[oaicite:0]{index=0}
Mesh.Recombine3DAll = 0;   // do NOT recombine into hex/prisms
Mesh.RecombineAll = 0;
Mesh.ElementOrder = 1;     // linear tets (Tet4)

Mesh 3;
Mesh.SaveAll = 1;
