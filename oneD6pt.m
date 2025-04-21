%% Matlab mesh
%% 1D, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 6;
msh.POS = [
0 0 0;
1 0 0;
0.1999999999995569 0 0;
0.3999999999989731 0 0;
0.599999999998945 0 0;
0.7999999999994725 0 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 1 3 4
 3 4 4
 4 5 4
 5 6 4
 6 2 4
];
msh.PNT =[
 1 2
 2 3
];
