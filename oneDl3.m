%% Matlab mesh
%% 1D, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 5;
msh.POS = [
0 0 0;
5 0 0;
1.249999999996376 0 0;
2.4999999999929 0 0;
3.749999999996682 0 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 1 3 4
 3 4 4
 4 5 4
 5 2 4
];
msh.PNT =[
 1 2
 2 3
];
