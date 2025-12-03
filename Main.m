%% Set-up inputs
clc,clear,close
addpath('functions');

%% INPUT:
%---------- WING GEOMETRY: ----------
c = 3;       %chord
b_span = 18; %wing span
N = 20;      %Number panels in chord direction
M = 20;      %Number panels in half wing (wing span direction) EVEN!!!
Sweep = deg2rad(0); %Sweep angle


%---------- FLOW: ----------
AoA = deg2rad(1);   %angle of attack
beta = deg2rad(0);  %angle of sideslip
U = 1;              %freestream intensity
rho = 1;            %density

%% Build geometry:
% Create U_infinity as vector
U_infinity = [U*cos(AoA)*cos(beta), -U*sin(beta), U*sin(AoA)*cos(beta)];

%------- WING: -------
geom = geometry(c,b_span,N,M,Sweep);

% Plot geometry 3D
plot3Dgeometry(geom)


%% SOLVE GAMMA:
[A,b] = scratc_system(geom,U_infinity,1);
GAMMA = A\b;
GAMMA = reshape(GAMMA,[2*M,N])';

%% POST PROCESSING:
%Compute all the aerdynamic loads:
[F,M,C_L,C_D,C_M,Cp] = aerodynamic_paramiters(geom,N,M,GAMMA,rho,U_infinity);

L_total = F(3);
D_total = F(1);
C_L_total = L_total/(0.5*rho*c*b_span*norm(U_infinity)^2);
C_D_total = D_total/(0.5*rho*c*b_span*norm(U_infinity)^2);
C_M_total = - M(2)/(0.5*rho*c^2*b_span*norm(U_infinity)^2);


%% Figures 
plotData = struct( ...
    'Cp', Cp, ...
    'gamma', GAMMA, ...
    'C_D', C_D, ...
    'C_L', C_L ...
);

plotFigures(plotData, geom)


%% VALIDATION (uncommend when you have XFLR5 data):
% Gamma_data = load('GG.txt');
% 
% figure()
% hold on
% plot(geom.centroids.y(1,:),GAMMA(1,:),'-b')
% plot(Gamma_data(:,1),Gamma_data(:,2),'kd', 'MarkerFaceColor', 'g')
% grid on
% legend('Code','Reference','fontsize',10)
% xlabel('Span','FontSize', 10,'fontweight','bold')
% ylabel('Gamma','FontSize', 10,'fontweight','bold')
% title('Gamma distribution along the span','FontSize', 15)
% hold off
% % saveas(gcf, 'GAMMA wing validation','png')
% 
% 
% 
% Cp_data = load('Cp_wing.txt');
% figure()
% hold on
% plot(geom.centroids.y(1,:),Cp(1,:),'-k')
% plot(geom.centroids.y(round(N/4),:),Cp(round(N/4),:),'-b')
% plot(geom.centroids.y(round(3*N/4),:),Cp(round(3*N/4),:),'-r')
% plot(Cp_data(27:39,1),Cp_data(27:39,2),'kd', 'MarkerFaceColor', 'g')
% plot(Cp_data(1:13,1),Cp_data(1:13,2),'kd', 'MarkerFaceColor', 'g')
% plot(Cp_data(14:26,1),Cp_data(14:26,2),'kd', 'MarkerFaceColor', 'g')
% grid on
% lg = legend('Leading edge','1/4 chord','3/4 chord','Reference','fontsize',10);
% set(lg,'location','best')
% xlabel('Span','FontSize', 10,'fontweight','bold')
% ylabel('C_p','FontSize', 10,'fontweight','bold')
% title('Pressure distribution along the span','FontSize', 15)
% hold off
% %saveas(gcf, 'Cp wing validation','png')
% 
% Cl_data = load('CL.txt');
% figure()
% hold on
% plot(geom.centroids.y(1,:),C_L(5,:),'-k')
% plot(Cl_data(:,1),Cl_data(:,2),'kd', 'MarkerFaceColor', 'g')
% grid on
% lg = legend('Code','Reference','fontsize',10);
% set(lg,'location','best')
% xlabel('Span','FontSize', 10,'fontweight','bold')
% ylabel('C_L','FontSize', 10,'fontweight','bold')
% title('Lift coefficient distribution along the span','FontSize', 15)
% hold off
% % saveas(gcf, 'Cp wing validation','png')
