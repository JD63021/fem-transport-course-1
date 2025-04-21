%% Matlab mesh
%% 1D, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 3;
msh.POS = [
0 0 0;
5 0 0;
2.4999999999929 0 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 1 3 4
 3 2 4
];
msh.PNT =[
 1 2
 2 3
];
