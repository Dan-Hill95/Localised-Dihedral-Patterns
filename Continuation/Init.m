global home
home=pwd; 
fprintf('%s\n',['Setting library path beginning with ' home]); 
addpath(genpath([home,'/Continuation_Codes']));