function [x, y, z, x_v, y_v, z_v, x_c, y_c, z_c, n,X_c,Y_c,Z_c] = geometry(c, b, N, M, sweepAngle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE THE GEOMETRY:

%INPUTS: chord (c), span (b), Number of panels in chord-direction (N) and
%half of the ones in the span-direction (M), sweep angle (sweepAngle)

%OUTPUTS: coordinates for extrema of the panels (x,y,z), extrema of horseshoe
%vortices (x_v,y_v,z_v) at 1/4 of the panel, control points (x_c,y_c,z_c)
%at 3/4 of the panel, the normal to each panel (n) and the coordinates of
%the center of each panel (X_c,Y_c,Z_c)

%Note: all the coordinates are given in a matricial form as if they
%rapresent the panel iteself (i.e. considering x(i,j),y(i,j),z(i,j) you have the
%point in position (i,j) of the panelled wing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Extrema of the panels:
x = linspace(0, c, N+1);
y = linspace(-b/2, b/2, 2*M+1);
z = zeros(1, (2 * M + 1) * (N + 1));
[x,y] = meshgrid(x,y);
x = x';
y = y';
z = reshape(z, [N+1, 2*M+1]);



% Panel length
k=1;
for i = 1:N
    for j = 1:2*M+1
        x_panel(k) = abs(x(i,j)-x(i+1,j)); % Length of panel in x-direction
        k = k+1;
    end
end
k=1;
for i =1:N+1
    for j = 1:2*M
        y_panel(k) = abs(y(i,j)-y(i,j+1)); % Length of panel in y-direction
        k = k+1;
    end
end
x_panel = reshape(x_panel,[2 * M+1, N])';
y_panel = reshape(y_panel,[2 * M, N+1])';




% Vortex extrema:
for i = 1:N
    for j = 1:(2*M+1)
        x_v(i,j) = x(i,j) + x_panel(i,j)/4;
        y_v(i,j) = y(i,j);
        z_v(i,j) = z(i,j);
    end
end



% Centroid coordinates:
for i = 1:N
    for j = 1:2*M
        x_c(i,j) = x(i,j)+3*x_panel(i,j)/4;
        y_c(i,j) = y(i,j)+y_panel(i,j)/2;
        z_c(i,j) = z(i,j);
    end
end



%Centre of panels:
for i = 1:N
    for j = 1:2*M
        X_c(i,j) = x(i,j)+x_panel(i,j)/2;
        Y_c(i,j) = y(i,j)+y_panel(i,j)/2;
        Z_c(i,j) = z(i,j);
    end
end



% Apply sweep angle to the wing
x(:,1:(size(x,2)/2)) = x(:,1:(size(x,2)/2)) - tan(sweepAngle) * y(:,1:(size(x,2)/2));
x(:,(size(x,2)/2+1):end) = x(:,(size(x,2)/2+1):end) + tan(sweepAngle) * y(:,(size(x,2)/2+1):end);
x_c(:,1:(size(x_c,2)/2)) = x_c(:,1:(size(x_c,2)/2)) - tan(sweepAngle) * y_c(:,1:(size(x_c,2)/2));
x_c(:,(size(x_c,2)/2+1):end) = x_c(:,(size(x_c,2)/2+1):end) + tan(sweepAngle) * y_c(:,(size(x_c,2)/2+1):end);
x_v(:,1:(size(x_v,2)/2)) = x_v(:,1:(size(x_v,2)/2)) - tan(sweepAngle) * y_v(:,1:(size(x_v,2)/2));
x_v(:,(size(x_v,2)/2+1):end) = x_v(:,(size(x_v,2)/2+1):end) + tan(sweepAngle) * y_v(:,(size(x_v,2)/2+1):end);
X_c(:,1:(size(X_c,2)/2)) = X_c(:,1:(size(X_c,2)/2)) - tan(sweepAngle) * Y_c(:,1:(size(X_c,2)/2));
X_c(:,(size(X_c,2)/2+1):end) = X_c(:,(size(X_c,2)/2+1):end) + tan(sweepAngle) * Y_c(:,(size(X_c,2)/2+1):end);


% Normal
n = [0, 0, 1]; %if you want to add a dihedral/twist angle you have to account for a different normal

end
