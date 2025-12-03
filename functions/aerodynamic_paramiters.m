function [F_tot,M_tot,C_L,C_D,C_M,Cp] = aerodynamic_paramiters(geom,N, M,GAMMA,rho,U_infinity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTE THE AERODYNAMICS LOADS:

%INPUTS: coordinates for extrema of the panels (x,y,z), extrema of horseshoe
%vortices (geom.vertices.x,geom.vertices.y,geom.vertices.z) at 1/4 of the panel, number of panels in the chord
%and span direction (M,N), circulation of each horseshoe vortex (GAMMA),
%air density (rho), freestream vector (U_infinity), coordinates of the
%center point of each panel (geom.centers.x,geom.centers.y,geom.centers.z)

%OUTPUTS:  Total Aerodynamics force and moment vector generated (F,Moment), 
%Lift Drag and Moment coefficient produced by each panel (C_L, C_D, C_M) and 
%Pressure coefficient of each panel (Cp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract geometric data
px = geom.panels.x;
py = geom.panels.y;
pz = geom.panels.z;

vx = geom.vortices.x;
vy = geom.vortices.y;
vz = geom.vortices.z;

cx = geom.centers.x;
cy = geom.centers.y;
cz = geom.centers.z;

%Define length of each panel in x and y direction
x_panel = abs(px(1:N, :) - px(2:N+1, :));    % (N) × (2M+1)
y_panel = abs(py(:, 1:2*M) - py(:, 2:2*M+1)); % (N+1) × (2M)

% Extrac constants
U = U_infinity;
U_mag2 = norm(U)^2;
rho2 = 0.5 * rho;

% Initialize vectors
F = zeros(N, 2*M, 3);
Moment = zeros(N, 2*M, 3);
C_L = zeros(N, 2*M);
C_D = zeros(N, 2*M);
C_M = zeros(N, 2*M);
Cp  = zeros(N, 2*M);
Lift = zeros(N, 2*M);
Drag = zeros(N, 2*M);

for i = 1:N
    for j = 1:2*M

        % Panel vortex endpoints
        v1 = [vx(i,j),   vy(i,j),   vz(i,j)];
        v2 = [vx(i,j+1), vy(i,j+1), vz(i,j+1)];

        % Panel center
        center = [cx(i,j), cy(i,j), cz(i,j)];

        % Induced velocity
        u_i = biotSavart(v1, v2, center, GAMMA(i,j));

        % Aerodynamic force
        F_ij = rho * GAMMA(i,j) * cross((U + u_i), [0, y_panel(i,j), 0]);
        F(i,j,:) = F_ij;

        % Force components
        Lift(i,j) = F_ij(3);
        Drag(i,j) = F_ij(1);

        % Coefficients
        A = x_panel(i,j) * y_panel(i,j);     % panel area
        C_L(i,j) = Lift(i,j) / (rho2 * A * U_mag2);
        C_D(i,j) = Drag(i,j) / (rho2 * A * U_mag2);

        % Moment about x-axis
        M_ij = cross(F_ij, [cx(i,j), 0, 0]);
        Moment(i,j,:) = M_ij;

        C_M(i,j) = -M_ij(2) / (rho2 * x_panel(i,j)^2 * y_panel(i,j) * U_mag2);
    end

%Total Aerodynamic force
F_tot = squeeze(sum(sum(F,1),2));
%Total Aerodynamic moment
M_tot = squeeze(sum(sum(Moment,1),2));

Cp = 2 * GAMMA ./ norm(U);

end
