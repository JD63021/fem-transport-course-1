// Gmsh project created on Mon Mar 24 12:14:19 2025
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {Pi, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Physical Point(2) = {1};
//+
Physical Point(3) = {2};
//+
Physical Curve(4) = {1};
