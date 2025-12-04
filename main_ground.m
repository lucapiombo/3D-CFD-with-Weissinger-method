%% Set-up inputs
clc,clear,close
addpath('functions');

%---------- WING GEOMETRY: ----------
c_wing = 3; %chord
b_wing = 21; %wing span
N_wing = 20; %Number panels in chord direction
M_wing = 20; %Number panels in half wing (wing span direction) EVEN!!!
Sweep_wing = deg2rad(0); %Sweep angle
h_wing = 3; %height

%---------- TAIL GEOMETRY: ----------
c_tail = 1.5; %chord
b_tail = 6; %wing span
N_tail = 20; %Number panels in chord direction
M_tail = 20; %Number panels in half wing (wing span direction) EVEN!!!
Sweep_tail = deg2rad(0); %Sweep angle
h_tail = 2.5;

%---------- DISTANCE TAIL-WING: ----------
x_wt = 5; %x-distance between wing-tail;
y_wt = 0; %distance LE_wing-LE_tail;

%---------- FLOW: ----------
alpha = deg2rad(1); %angle of attack
beta = deg2rad(0); %angle of sideslip
U = 1; %freestream intensity
rho = 1.225; %density

%% Build geometry:

%Create U_infinity as vector
U_infinity = [U*cos(alpha)*cos(beta), -U*sin(beta), U*sin(alpha)*cos(beta)]; %U_infinity vector

% Create geometry
[wing, imgWing] = geometry(c_wing,b_wing,N_wing,M_wing,Sweep_wing, z=h_wing);
[tail, imgTail] = geometry(c_tail,b_tail,N_tail,M_tail,Sweep_tail, z=h_tail, x=x_wt, y=y_wt);

% Plot geometry 3D
plot3Dgeometry({wing, imgWing, tail, imgTail})

%% SOLVE GAMMA:

[A_ww,b_w] = scratc_system(wing, wing, U_infinity, 1);
[A_ww_im] = scratc_system(wing, imgWing, U_infinity, 1);
A_ww = A_ww-A_ww_im;

[A_tw] = scratc_system(wing, tail, U_infinity, 1);
[A_tw_im] = scratc_system(wing, imgTail, U_infinity, 1);
A_tw = A_tw+A_tw_im;

[A_tt,b_t] = scratc_system(tail, tail, U_infinity, 1);
[A_tt_im] = scratc_system(tail, imgTail, U_infinity, 1);
A_tt = A_tt+A_tt_im;

[A_wt] = scratc_system(tail, wing, U_infinity, 1);
[A_wt_im] = scratc_system(tail, imgWing, U_infinity, 1);
A_wt = A_wt+A_wt_im;

% Build up the system:
A = zeros(2*M_wing*N_wing+2*M_tail*N_tail);
A(1:(2*M_wing*N_wing),1:(2*M_wing*N_wing)) = A_ww;
A(1:(2*M_wing*N_wing),(2*M_wing*N_wing+1):end) = A_tw;
A((2*M_wing*N_wing+1):end,1:(2*M_wing*N_wing)) = A_wt;
A((2*M_wing*N_wing+1):end,(2*M_wing*N_wing+1):end) = A_tt;
b = [b_w; b_t];

GAMMA = A\b;

GAMMA_w = reshape(GAMMA(1:(2*M_wing*N_wing))', [2*M_wing, N_wing])';
GAMMA_t = reshape(GAMMA((2*M_wing*N_wing+1):end)', [2*M_tail, N_tail])';

%% POST PROCESSING:

[F_w,Moment_w,C_Lw,C_Dw,C_Mw,Cp_w] = aerodynamic_paramiters(wing, N_wing, ...
                                                            M_wing, GAMMA_w, ...
                                                            rho, U_infinity);
[F_t,Moment_t,C_Lt,C_Dt,C_Mt,Cp_t] = aerodynamic_paramiters(tail, N_tail, ...
                                                            M_tail, GAMMA_t, ...
                                                            rho, U_infinity);

L_total_w = F_w(3);
D_total_w = F_w(1);
C_L_total_w = L_total_w/(0.5*rho*c_wing*b_wing*norm(U_infinity)^2);
C_D_total_w = D_total_w/(0.5*rho*c_wing*b_wing*norm(U_infinity)^2);
C_M_total_w = - Moment_w(2)/(0.5*rho*c_wing^2*b_wing*norm(U_infinity)^2);

L_total_t = F_t(3);
D_total_t = F_t(1);
C_L_total_t = L_total_t/(0.5*rho*c_tail*b_tail*norm(U_infinity)^2);
C_D_total_t = D_total_t/(0.5*rho*c_tail*b_tail*norm(U_infinity)^2);
C_M_total_t = - Moment_t(2)/(0.5*rho*c_tail^2*b_tail*norm(U_infinity)^2);

%% FIGURES:
plotData_w = struct( ...
    'Cp', Cp_w, ...
    'gamma', GAMMA_w, ...
    'C_D', C_Dw, ...
    'C_L', C_Lw ...
);

plotData_t = struct( ...
    'Cp', Cp_t, ...
    'gamma', GAMMA_t, ...
    'C_D', C_Dt, ...
    'C_L', C_Lt ...
);

plotFigures({plotData_w, plotData_t}, {wing, tail})