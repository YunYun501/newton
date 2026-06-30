SetFactory("OpenCASCADE");

// Coarse version (lc=6mm) for Warp FEM solver
R  = 19;
L  = 49;
t  = 2.5;
lc = 6.0;
eps = 1e-3;

Mesh.CharacteristicLengthMin = lc * 0.5;
Mesh.CharacteristicLengthMax = lc * 1.5;

Cylinder(1) = {0, 0, 0,  0, 0, L,  R};
Box(2) = {-eps, -R-eps, -eps,  R+2*eps, 2*(R+eps), L+2*eps};
outerHalf[] = BooleanIntersection{ Volume{1}; Delete; }{ Volume{2}; Delete; };

Rin = R - t;
Lin = L;
Cylinder(3) = {0, 0, 0,  0, 0, Lin,  Rin};
Box(4) = {-eps, -Rin-eps, -eps,  Rin+2*eps, 2*(Rin+eps), Lin+2*eps};
innerHalf[] = BooleanIntersection{ Volume{3}; Delete; }{ Volume{4}; Delete; };

solid[] = BooleanDifference{ Volume{outerHalf[]}; Delete; }{ Volume{innerHalf[]}; Delete; };

Physical Volume("Chamber") = {solid[]};

Mesh.Algorithm3D = 1;
Mesh.Recombine3DAll = 0;
Mesh.ElementOrder = 1;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Mesh 3;
Mesh.SaveAll = 1;
