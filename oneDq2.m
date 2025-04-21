%% Matlab mesh
%% 1D, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 9;
msh.POS = [
0 0 0;
5 0 0;
1.249999999996376 0 0;
2.4999999999929 0 0;
3.749999999996682 0 0;
0.6249999999984707 0 0;
1.874999999994393 0 0;
3.124999999994791 0 0;
4.374999999998949 0 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES3 =[
 1 3 6 4
 3 4 7 4
 4 5 8 4
 5 2 9 4
];
msh.PNT =[
 1 2
 2 3
];
