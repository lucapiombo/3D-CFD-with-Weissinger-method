%% Set-up inputs
clc,clear,close
addpath('functions');

%---------- WING GEOMETRY: ----------
c = 3;       %chord
b_span = 18; %wing span
N = 5;      %Number panels in chord direction
M = 6;      %Number panels in half wing (wing span direction) EVEN!!!
Sweep = deg2rad(0); %Sweep angle

%---------- FLOW: ----------
AoA = deg2rad(1);   %angle of attack
beta = deg2rad(0);  %angle of sideslip
U = 1;              %freestream intensity
rho = 1;            %density

%% Build geometry:
% Create U_infinity as vector
U_infinity = [U*cos(AoA)*cos(beta), -U*sin(beta), U*sin(AoA)*cos(beta)];

% Create geometry
wing = geometry(c,b_span,N,M,Sweep);

% Plot geometry 3D
plot3Dgeometry(wing)

%% SOLVE GAMMA:
[A,b] = scratc_system(wing,U_infinity,1);
GAMMA = A\b;
GAMMA = reshape(GAMMA,[2*M,N])';

%% POST PROCESSING:
[F,M,C_L,C_D,C_M,Cp] = aerodynamic_paramiters(wing,N,M,GAMMA,rho,U_infinity);

L_total = F(3);
D_total = F(1);
C_L_total = L_total/(0.5*rho*c*b_span*norm(U_infinity)^2);
C_D_total = D_total/(0.5*rho*c*b_span*norm(U_infinity)^2);
C_M_total = - M(2)/(0.5*rho*c^2*b_span*norm(U_infinity)^2);

%% FIGURES:
plotData = struct( ...
    'Cp', Cp, ...
    'gamma', GAMMA, ...
    'C_D', C_D, ...
    'C_L', C_L ...
);

plotFigures(plotData, wing)