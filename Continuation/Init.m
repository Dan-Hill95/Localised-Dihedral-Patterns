global home
home=pwd; 
fprintf('%s\n',['Setting library path beginning with ' home]); 
addpath(genpath([home,'/Continuation_Codes']));
fprintf('%s\n',[' ']); 
fprintf('%s\n',['These codes are for the continuation of localised dihedral patterns in the 2-3 Swift-Hohenberg equation.']); 
fprintf('%s\n',['Here are your next steps:']); 
fprintf('%s\n',['1. If you want to solve the (n+1)-dimensional algebraic matching condition for localised dihedral']); 
fprintf('%s\n',['   patterns, type:']); 
fprintf('%s\n',[' ']); 
fprintf('%s\n',['>> a = MatchSoln(x, m, r_max, mu);']); 
fprintf('%s\n',[' ']); 
fprintf('%s\n',['   where x=[x[0],x[1],...,x[n]] is your initial guess, m is the dihedral lattice (must be even!),']); 
fprintf('%s\n',['   r_max is the radial boundary and mu gives an approximate localisation. For example:']); 
fprintf('%s\n',[' ']); 
fprintf('%s\n',['>> a = MatchSoln((1/10)*ones(1,11), 6, 100, 1e-4);']); 
fprintf('%s\n',[' ']); 
fprintf('%s\n',['2. If you want to find and continue localised dihedral patterns in the 2-3 SH equation, type:']); 
fprintf('%s\n',[' ']); 
fprintf('%s\n',['>> branch = Cont_Patch(x, p, Dir);']); 
fprintf('%s\n',[' ']); 
fprintf('%s\n',['   where p=[mu, gamma, kappa, m, n+1] contains the parameters of the system including the bifurcation']); 
fprintf('%s\n',['   parameter mu, the respective quadratic and cubic coefficients (gamma, kappa), and the dimension ']); 
fprintf('%s\n',['   (n+1) of the reduced ODE system. Dir indicates the direction and step size of the continuation']); 
fprintf('%s\n',['   routine, which must be ''pl'' (plus), ''mn'' (minus), ''sp'' (small plus), or ''sm'' (small minus).']); 
fprintf('%s\n',['   For example:']); 
fprintf('%s\n',[' ']); 
fprintf('%s\n',['>> branch = Cont_Patch(1*[-1,2], [0.02, 1.6, -1, 2, 5], ''pl'');']); 
fprintf('%s\n',[' ']); 
fprintf('%s\n',['3. If you want to plot the bifurcation curve of a localised solution, type:']); 
fprintf('%s\n',[' ']); 
fprintf('%s\n',['>> ExploreBifurcationDiagram(''FolderName/branch.mat'',idVar);']); 
fprintf('%s\n',[' ']); 
fprintf('%s\n',['   where FolderName is the data folder that is created by Cont_Patch, and idVar is the column of the']); 
fprintf('%s\n',['   branch measures to be plotted. For example, for a data folder called ''D2_Patch_pl'', you would type:']); 
fprintf('%s\n',[' ']); 
fprintf('%s\n',['>> ExploreBifurcationDiagram(''D2_Patch_pl/branch.mat'',5);']); 
fprintf('%s\n',[' ']); 