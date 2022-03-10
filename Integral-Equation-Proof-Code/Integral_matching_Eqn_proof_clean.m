% Rigorous proof for a solution to the integral equation
% w(t) = 2\int_t^1 w(s)w(s-t)ds + int_0^t w(s)w(t-s)ds
% the code accompanies the paper
% Approximate Localised Dihedral Patterns Near a Turing Instability
% by DJ Hill, JJ Bramburger, DJB Lloyd
%
% code follows that of 
% Rigorous computation of a radially symmetric localised solution in a Ginzburg-Landau problem
% by J.B van den Berg, C.M. Groothedde and J.F. Williams, SIADS 
% with various modifications
%
% The approximate solution to G^m(w)=0 and the chosen mesh are stored in 
% the Integral_proof_N1001_data.mat.
% This approximate solution was obtained by using Newton's method on the
% trapezoidal approximation of the integral equation
% 
% Before running one needs to initialise IntLab:
% https://www.tuhh.de/ti3/rump/intlab/

N = 1001; % mesh size for initial guess
load Integral_proof_N1001_data.mat;

%startintlab; % run intlab
t = intval(data(:,1));
wh = intval(data(:,2));
n = length(wh)-1; %Number of intervals
NN= intval(N);
omega = intval('0.02'); 

fprintf('Using %i mesh points.\n',length(t));
fprintf('Using omega=%d.\n',mid(omega));

w1 = intval(wh(1:n)); %Endpoints of the intervals [w(t_k), w(t_(k+1))]
w2 = intval(wh(2:n+1));
t1 = intval(t(1:n)); %Endpoints of the intervals [w(t_k), w(t_(k+1))]
t2 = intval(t(2:n+1));
dt = intval(sup(max(t2-t1)));

%% calculate convolution-type integrals
F=zeros(N,1); S1=zeros(N,1); S2=zeros(N,1);

%Introduce correlation and convolution operators
T1 = intval(toeplitz([flipud(wh) ; zeros((N-1),1)],zeros(N,1)));
T2 = intval(toeplitz([wh ; zeros((N-1),1)],zeros(N,1)));
temp1 = T1*wh; temp2 = T2*wh;

%Extract the relevant part of each operator
S1=temp1(N:2*N-1); S2=temp2(1:N);

%Define Gm
Gm = intval((1+2*wh(1)/NN)*wh - 2*S1/NN - S2/NN + wh(N)*flipud(wh)/NN);

fprintf('Gm is calculated.\n');

% DGm
e=sparse(N,N); e(1,1)=1; e = intval(e); %Matrix version of unit vector in a_0 direction
J0 = intval(flipud(speye(N))); %Define exchange matrix (sparse matrix with ones on main anti-diagonal)
Temp1=T2*J0 + T1; %Linearisation of first sum
Temp2=2*T2; %Linearisation of second sum
dS1=Temp1(N:2*N-1,1:N);
dS2=Temp2(1:N,1:N); 

DGm = intval((1+2*wh(1)/NN)*eye(N) + 2*diag(wh)*e/NN - wh(N)*eye(N)/NN - diag(flipud(wh))*flipud(e)/NN - 2*dS1/NN - dS2/NN);   %Define the Jacobian
fprintf('DGm is calculated.\n');
 
%Next we compute the actual radii-polynomials. We first calculate the
%Y-bounds, followed by the Z bounds.
Ah = intval(inv(mid(DGm))); % approximate 
fprintf('Ah is calculated.\n');
 
Y_k = abs(Ah*Gm);
fprintf('Y_k is calculated.\n');

%We define for each k intervals with possible ranges for s and the 
%piece-wise linear function wh. 
s_I = infsup(inf(t(1:n)),sup(t(2:n+1)));
wh_I = infsup(min([inf(wh(2:n+1)) inf(wh(1:n))]')',max([sup(wh(2:n+1)) sup(wh(1:n))]')'); 
dwh_I = infsup(min([inf(wh(2:n+1))-inf(wh(1:n)) inf(wh(2:n+1))-sup(wh(1:n)) sup(wh(2:n+1))-inf(wh(1:n)) sup(wh(2:n+1))-sup(wh(1:n)) ]')',...
               max([inf(wh(2:n+1))-inf(wh(1:n)) inf(wh(2:n+1))-sup(wh(1:n)) sup(wh(2:n+1))-inf(wh(1:n)) sup(wh(2:n+1))-sup(wh(1:n)) ]')');

% Y infinity bound
h3 = 2*wh_I(end:-1:1).*dwh_I(end) + 2*dwh_I(end:-1:1).*wh_I(n) ...
    + dwh_I*wh_I(1) + wh_I*dwh_I(1);

Y_inf = max(sup(((diam(s_I))).*(abs(h3))/8));
fprintf('Y_inf is calculated.\n');

% Z_k bounds
V1 = intval(zeros(N,1));
V1(1) = dt*([1 2*ones(1,N-2) 1]*wh); % V1_{k=0}
V1(2) = dt/2*([1 2*ones(1,N-3) 1 0]*wh + (wh(1) + wh(2)) ) + V1(1)/2; % V1_{k=1}
for i = 3:N-2
    V1(i) = dt/2*([1 2*ones(1,N-i-1) 1 zeros(1,i-1)]*wh ...
            + [1 2*ones(1,i-2) 1 zeros(1,N-i)]*wh) + V1(1)/2;
end
V1(N-1) = V1(2); V1(N) = V1(1);

V2 = intval(zeros(n+1,1));
V2(:) = abs(intval(2*(1+omega)^2));

Z_k = intval(zeros(n+1,2));
Z_k(:,1) = (abs(eye(n+1)-Ah*DGm))*ones(n+1,1) +  abs(Ah)*V1*omega; 
Z_k(:,2) = abs(Ah)*V2;
fprintf('Z_k is calculated.\n');

% Z_infty bounds
W1 = max(sup(4*(1+omega)*diam(s_I)*(2*sup(wh_I(1)) + sup(wh_I(end)))));
W2 = max(sup(6*(1+omega)*(2*diam(s_I)*(1+omega)+omega)));

Z_inf = intval(zeros(1,2));
Z_inf(1) = W1;
Z_inf(2) = W2; 
fprintf('Z_inf is calculated.\nNow calculating and solving radii-polynomials.\n');

%The following are the actual coefficients of the radii-polynomials
rad_p = [Z_k(:,2) Z_k(:,1)-ones(length(Z_k),1) Y_k; Z_inf(:,2) Z_inf(:,1)-omega Y_inf];
 
%The following part calculates the parts where polynomials are negative.
rad_I = zeros(n+1,2);
for k=1:n+2
    rts=sort(roots(sup(rad_p(k,:))));
    rad_I(k,:) = rts(1:2);
end
 
min_r = min(rad_I) %zeros must lie below the min of the upper bounds
max_r = max(rad_I) %zeros must lie above the max of the lower bounds
 
fprintf('Taking r=0.05...\n');
r = 0.05;
if max(sup((rad_p*[r^2; r; 1])))<0
    fprintf('Success! Largest value of p(r) is %d\nAll radii-polynomials are negative!\n', max(sup((rad_p*[r^2; r; 1]))));
end
if min(inf(wh))>sup(r*(1+omega))
    fprintf('Furthermore, min(wh)>r(1+omega), hence the solution is positive.\n')
end
