LEN = 0.005;
WID = 0.0025;
THI = 0.0005;
RAD = 0.0002;
SPC = 0.0002;

Y_SPACE = WID / 5;

Point(1) = {0.0, 0.0, 0.0};
Extrude {LEN/2, 0.0, 0.0} { Point{1}; }
Extrude {0.0, WID/2, 0.0} { Line{1}; }
Delete { Surface{5}; }
Extrude {0.0, 0.0, THI/2} { Line{1:4}; }

Point(10) =          {SPC+2*RAD, WID/2-1*Y_SPACE, 0.0};
Extrude { {0, 0, 1}, {SPC+1*RAD, WID/2-1*Y_SPACE, 0.0}, Pi/2} { Point{10}; }
Extrude { {0, 0, 1}, {SPC+1*RAD, WID/2-1*Y_SPACE, 0.0}, Pi/2} { Point{11}; }
Extrude { {0, 0, 1}, {SPC+1*RAD, WID/2-1*Y_SPACE, 0.0}, Pi/2} { Point{13}; }
Extrude { {0, 0, 1}, {SPC+1*RAD, WID/2-1*Y_SPACE, 0.0}, Pi/2} { Point{14}; }
Extrude {0.0, 0.0, THI/2} { Curve{21:24}; }

Point(30) =          {SPC+2*RAD, WID/2-2*Y_SPACE, 0.0};
Extrude { {0, 0, 1}, {SPC+1*RAD, WID/2-2*Y_SPACE, 0.0}, Pi/2} { Point{30}; }
Extrude { {0, 0, 1}, {SPC+1*RAD, WID/2-2*Y_SPACE, 0.0}, Pi/2} { Point{31}; }
Extrude { {0, 0, 1}, {SPC+1*RAD, WID/2-2*Y_SPACE, 0.0}, Pi/2} { Point{33}; }
Extrude { {0, 0, 1}, {SPC+1*RAD, WID/2-2*Y_SPACE, 0.0}, Pi/2} { Point{34}; }
Extrude {0.0, 0.0, THI/2} { Curve{41:44}; }

Curve Loop(1) = {1, 4, -2, -3};
Curve Loop(2) = {44, 41, 42, 43};
Curve Loop(3) = {24, 21, 22, 23};
Plane Surface(61) = {1, 2, 3};

Curve Loop(4) = {5, 17, -9, -13};
Curve Loop(5) = {53, 57, 45, 49};
Curve Loop(6) = {37, 25, 29, 33};
Plane Surface(62) = {4, 5, 6};

Surface Loop(1) = {62, 8, 61, 20, 12, 16, 56, 60, 48, 52, 36, 40, 28, 32};
Volume(1) = {1};

Physical Volume(0) = {1};

Physical Surface(0) = {20}; // prescribe zero displacement on X axis
Physical Surface(1) = {8};  // prescribe zero displacement on Y axis
Physical Surface(2) = {61}; // prescribe zero displacement on Z axis
Physical Surface(6) = {12,16,28,40,48,60,62}; // traction free
Physical Surface(10) = {32,36,52,56}; // prescribe non-zero displacement on Z axis


Mesh.Algorithm3D = 4;
Mesh.OptimizeNetgen = 1;
Mesh.CharacteristicLengthMin = 0.00002;
Mesh.CharacteristicLengthMax = 0.00010;
Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;
