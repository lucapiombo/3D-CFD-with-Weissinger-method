function [F,Moment,C_L,C_D,C_M,Cp] = aerodynamic_paramiters(x,y,x_v,y_v,z_v,N, M,GAMMA,rho,U_infinity,X_c,Y_c,Z_c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTE THE AERODYNAMICS LOADS:

%INPUTS: coordinates for extrema of the panels (x,y,z), extrema of horseshoe
%vortices (x_v,y_v,z_v) at 1/4 of the panel, number of panels in the chord
%and span direction (M,N), circulation of each horseshoe vortex (GAMMA),
%air density (rho), freestream vector (U_infinity), coordinates of the
%center point of each panel (X_c,Y_c,Z_c)

%OUTPUTS:  Total Aerodynamics force and moment vector generated (F,Moment), 
%Lift Drag and Moment coefficient produced by each panel (C_L, C_D, C_M) and 
%Pressure coefficient of each panel (Cp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define length of each panel in x and y direction
k=1;
for i = 1:N
    for j = 1:2*M+1
        x_panel(k) = abs(x(i,j)-x(i+1,j));
        k = k+1;
    end
end
k=1;
for i =1:N+1
    for j = 1:2*M
        y_panel(k) = abs(y(i,j)-y(i,j+1));
        k = k+1;
    end
end
x_panel = reshape(x_panel,[2 * M+1, N])';
y_panel = reshape(y_panel,[2 * M, N+1])';


% Compute F, Moment, C_L, C_D, C_M:
k = 1;
F = zeros(2*M*N,3);
Moment = zeros(2*M*N,3);
for i = 1:N
    for j = 1:(2*M)
%         u_i = [0,0,0];
%         for l = 1:N
%             for m = 1:(2*M)
                u_i = Biot_Savart([x_v(i,j),y_v(i,j),z_v(i,j)],[x_v(i,j+1),y_v(i,j+1),z_v(i,j+1)],[X_c(i,j),Y_c(i,j),Z_c(i,j)],GAMMA(i,j));
%             end
%         end
        F(k,:) = rho*GAMMA(i,j)*cross((U_infinity+u_i), [0,y_panel(i,j),0]); %Aerodynamic force

        Lift(i,j) = F(k,3); %Lift
        C_L(i,j) = Lift(i,j)/(0.5*rho*x_panel(i,j)*y_panel(i,j)*norm(U_infinity)^2); %Lift coefficient
        
        Drag(i,j) = F(k,1); %Drag
        C_D(i,j) = Drag(i,j)/(0.5*rho*x_panel(i,j)*y_panel(i,j)*norm(U_infinity)^2); %Drag coefficient

        Moment(k,:) = cross(F(k,:),[X_c(i,j),0,0]); %Aerodynamic moment
        C_M(i,j) = -Moment(k,2)/(0.5*rho*x_panel(i,j)^2*y_panel(i,j)*norm(U_infinity)^2); %Moment coefficient

        k = k+1;
    end
end
F = sum(F); %Total Aerodynamic force
Moment = sum(Moment); %Total Aerodynamic moment


%Cp:
for i = 1:N
    for j = 1:2*M
        Cp(i,j) = 2*GAMMA(i,j)/norm(U_infinity);
    end
end

% %Check if the central points are correct:
% figure
% surf(x, y, z)
% hold on
% plot3(X_c, Y_c, Z_c,'*b')
% hold off
end
