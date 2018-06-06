//SetFactory("OpenCASCADE");

a = 0.1;
clcrack = a/5.0;
// By changing this we can make it nonuniform
// We keep it uniform to be fair to peridynamics, which doesn't
// work on non uniform grid!
//clbig = clcrack;
clbig = 0.5;

// Rectangle(1) = {-1, -1, 0, 1, 1, 0};
Point(1) = {-1,-1,0,clbig};
Point(2) = { 1,-1,0,clbig};
Point(3) = { 1, 1,0,clbig};
Point(4) = {-1, 1,0,clbig};

Point(5) = {-1, 0,0,clbig};
Point(6) = {-a, 0,0,clcrack};
Point(7) = { a, 0,0,clcrack};
Point(8) = { 1, 0,0,clbig};


//  4<-------3
//  v        ^
//  5->6->7->8
//  v        ^
//  1------->2

Line(1) = {1,2};
Line(2) = {2,8};
Line(3) = {8,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,1};
Line(7) = {5,6};
Line(8) = {6,7}; // the crack
Line(9) = {7,8};

Line Loop(1) = {1,2, -9,-8,-7,6};
Line Loop(2) = {7,8,9,3,4,5};

Plane Surface(1) = { 1 };
Plane Surface(2) = { 2 };

Physical Surface(1) = {1,2};
Physical Line(2) = { 8 };
Physical Line(3) = { 1,2,6,3,4,5 };

Mesh 2;

Plugin(Crack).Dimension = 1;
Plugin(Crack).PhysicalGroup = 2;
Plugin(Crack).OpenBoundaryPhysicalGroup = 0;
Plugin(Crack).Run ;

Physical Line(4) = { 10 };

Save "crack2.msh";
