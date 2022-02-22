function [runtime,w_old] = vonHardenbergLOD2D_IFETDRDP(dt,steps)

% dt: time step
% steps: number of spatial points in each coordinate direction

%% Model Paramters and initial conditions
G = 1.6; S = 1.6; R = 1.5; N = 0.2; B = 3;

pc = N/(G-S*N); 

% diffusion coefficient
epsln1 = 1; epsln2 = 27;
D = epsln2;


K = 0.234;
p = 0.223;
P = 0.04-(p-pc)/D;

% create nodes
x = linspace(-200,200,steps); 
h = abs(x(1)-x(2));
fprintf('dt=%f k*h=%f\n',dt,K*h)
y = x;
nnodes = steps^2;
nodes = zeros(nnodes,2);
j = 1;
for k = 1 : steps
        for i = 1:steps
               nodes(j,:) = [x(i) y(k)];
            j = j+1;
        end
end
nb = 2*nnodes; % becuase we are solving a syste of 2 RDE

% discretize time interval
t = 0:dt:500; tlen = length(t);

%Localised hexagon shape
hex=2*(cos(K*nodes(:,1))+ cos(K*(nodes(:,1)+sqrt(3)*nodes(:,2))/2)+ cos(K*(nodes(:,1)-sqrt(3)*nodes(:,2))/2))/3;

% initial condition for u
u_old = 0.05 + 0.3*exp(-sqrt(P)*(sqrt(nodes(:,1).^2+nodes(:,2).^2))/2).*hex; 

% initial condition for v
v_old =p + 0.05 + 0.3*exp(-sqrt(P)*(sqrt(nodes(:,1).^2+nodes(:,2).^2))/2).*hex; 


% Stacking nodes for evolution
w_old = zeros(nb,1);
w_old(1:2:nb-1) = u_old; w_old(2:2:nb) = v_old - (B*D/(D-1))*u_old; 

%% Block matrix Assembly
C = zeros(2);
C(1,1) = (epsln1*dt)/h^2;
C(2,2) = (epsln2*dt)/h^2;
Q = blktridiag(2,-1,-1,steps);
Q(1,2) = -2; Q(steps,steps-1) = -2; I = speye(steps);
A1 = kron(kron(I,Q),C);A2 = kron(kron(Q,I),C); 
I = speye(nb); 
M1 = sparse(I+A1); M2 = sparse(I+(1/3)*A1); M3 = sparse(I+(1/4)*A1);
M11 = sparse(I+A2); M22 = sparse(I+(1/3)*A2); M33 = sparse(I+(1/4)*A2);

%% Time Evolution 
[L1,U1] = lu(M1);
[L2,U2] = lu(M2);
[L3,U3] = lu(M3);
[L11,U11] = lu(M11); 
[L22,U22] = lu(M22);
[L33,U33] = lu(M33);

 Usoln = reshape(u_old,steps,steps); 
 Vsoln = reshape(v_old+(B*D/(D-1))*u_old,steps,steps);
 Uplot = Usoln(:,:);
 Vplot = Vsoln(:,:);
     scrsz = get(0,'ScreenSize');  
    
    figure(1);close gcf;figure(2);close gcf; 
    map=imadjust(flipud(summer),[.2 .3 0; .6 .7 1],[]); 
    
    
    figure('Position',[1*scrsz(3)/4 scrsz(4)/2 scrsz(3)/4 scrsz(4)/3]);
    surf(x,y,Uplot')
    xlabel('x')
    ylabel('y')
    colormap(map);colorbar;
    set(gca,'LineWidth', 1);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    pbaspect(gca,[1 1 1]);view(0,90);shading interp;
    set(gca,'XTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*(x(end)));
    set(gca,'YTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*(x(end)));
    
    figure('Position',[2*scrsz(3)/4 scrsz(4)/2 scrsz(3)/4 scrsz(4)/3]);
    surf(x,y,Vplot')
    xlabel('x')
    ylabel('y')
    colormap(flipud(parula));colorbar;
    set(gca,'LineWidth', 1);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    pbaspect(gca,[1 1 1]);view(0,90);shading interp;
    set(gca,'XTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*(x(end)));
    set(gca,'YTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*(x(end)));
 pause;
tic
for i = 2:tlen
     F_old = F(w_old);
     w_star = U11\(L11\(w_old + dt*F_old));
     w_star = U1\(L1\w_star);
     F_star = F(w_star);
     a_1 = U2\(L2\w_old);
     b_1 = U3\(L3\w_old);
     c_1 = 9*a_1 - 8*b_1;
     a_2 = U2\(L2\F_old);
     b_2 = U3\(L3\F_old);
     c_2 = 9*a_2-8*b_2;
     d_1 = U22\(L22\(9*c_1+2*dt*c_2 + dt*F_star));
     d_2 = U33\(L33\(-8*c_1-(3/2)*dt*c_2-(dt/2)*F_star));
     w_old = d_1+d_2;              
     
    if mod(i*dt,1)==0 
    u_soln = w_old(1:2:nb-1); 
    v_soln = w_old(2:2:nb)+(B*D/(D-1))*w_old(1:2:nb-1);
    Usoln = reshape(u_soln,steps,steps); 
    Vsoln = reshape(v_soln,steps,steps);
    Uplot = Usoln(:,:);
    Vplot = Vsoln(:,:);
        
    figure(1);
    surf(x,y,Uplot')
    xlabel('x')
    ylabel('y')
    colormap(map);colorbar;
    set(gca,'LineWidth', 1);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    pbaspect(gca,[1 1 1]);view(0,90);shading interp;
    set(gca,'XTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*(x(end)));
    set(gca,'YTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*(x(end)));
    
    figure(2);
    surf(x,y,Vplot')
    xlabel('x')
    ylabel('y')
    colormap(flipud(parula));colorbar;
    set(gca,'LineWidth', 1);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    pbaspect(gca,[1 1 1]);view(0,90);shading interp;
    set(gca,'XTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*(x(end)));
    set(gca,'YTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*(x(end)));
    end
end

runtime = toc;

function Fr = F(U)
 Fr = zeros(nb,1);
 u = U(1:2:nb-1); v = U(2:2:nb)+(B*D/(D-1))*U(1:2:nb-1);
%  f1 = A+u.^2.*v -(B+1)*u;
%  f2 = B*u-u.^2.*v;
 f1 = ((G*v)./(1+S*v)).*u - u.^2 - N*u;
 f2 = p - (1-R*u).*v-u.*v.^2;
 Fr(1:2:nb-1) = f1; Fr(2:2:nb) = f2 - (B*D/(D-1))*f1;
end

end
