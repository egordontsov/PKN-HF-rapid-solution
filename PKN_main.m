%Rapid approximate solution for PKN fracture
%https://arxiv.org/abs/2110.13088
%Egor Dontsov

clear all;clc;

%input parameters (do not make exactly zero to avoid division by zero)
E = 20;%GPa Young's modulus
nu = 0.2;%Poisson's ratio
KIc = 1;%MPa*m^1/2 fracture toughness
mu = 0.01;%Pa*s fluid viscosity
Cl = 1e-3;%mm/s^1/2 leakoff coefficient
H = 20;%m fracture height
Q0 = 1;%l/s injection rate
t = 1000;%s injection time

%set mesh promerties
Nt = 100;%number of time steps
Nx = 100;%number of spatial points

%scale problem parameters
Cp = 2*Cl;
Ep = E/(1-nu^2);

%run fast PKN solver
t = linspace(t/Nt,t,Nt)';
xi = linspace(0,1,Nx)';
[wvst,wvsx,lvst,etavst] = FastPKNSolver(t,xi,Cp,Ep,KIc,mu,H,Q0);

%load vertex solutions
[Wm,Lm,Wmt,Lmt,Wk,Lk,Wkt,Lkt] = PKNVertexSolutions(t,xi,Cp,Ep,KIc,mu,H,Q0);

%select vertex to plot
lvert = Lm;
wvert = Wm;
col = 'b';

%plot length versus time
figure;
plot([0;t],[0;lvst],'k-','linewidth',1.5);
hold on;
plot([0;t],[0;lvert],'--','color',col,'linewidth',2);
xlabel('t [s]','fontsize',16);
ylabel('l [m]','fontsize',16);

%plot efficency versus time
figure;
plot([0;t],[1;etavst],'k-','linewidth',1.5);
xlabel('t [s]','fontsize',16);
ylabel('\eta','fontsize',16);

%plot spatial width variation (this is height averaged width, multiply the result by 4/pi to get width at the center)
figure;
plot([xi*lvst(end); lvst(end)],[wvsx;0],'k-','linewidth',1.5);
hold on;
plot([xi*lvert(end); lvert(end)],[wvert;0],'--','color',col,'linewidth',2);
xlabel('x [m]','fontsize',16);
ylabel('w [mm]','fontsize',16);

%plot parametric space
PlotPKNParametricSpace(t,Cp,Ep,KIc,mu,H,Q0);
