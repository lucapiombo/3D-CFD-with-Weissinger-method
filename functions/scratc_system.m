function [A,b] = scratc_system(geom,U_infinity,GAMMA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTE THE SYSTEM COEFFICIENT:

%INPUTS: coordinates for extrema of the panels (x,y,z), extrema of horseshoe
%vortices (geom.vortices.x,geom.vortices.y,geom.vortices.z) at 1/4 of the panel, control points (geom.centroids.x,geom.centroids.y,geom.centroids.z)
%at 3/4 of the panel, the normal to each panel (n), streamflow vector 
%(U_infinity) and the gamma of the vortex (GAMMA)

%OUTPUTS: Coefficient of the A-matrix (A) and b-vector (b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract variables
n = geom.normal;
x_c = geom.centroids.x(:);
y_c = geom.centroids.y(:);
z_c = geom.centroids.z(:);

x_v = geom.vortices.x;
y_v = geom.vortices.y;
z_v = geom.vortices.z;


x_c = reshape(x_c',1,[]);
y_c = reshape(y_c',1,[]);
z_c = reshape(z_c',1,[]);


for i = 1:(size(x_c,2)*size(x_c,1))
    l = 1;
    k = 1;
    for j = 1:((size(x_v,2)-1)*size(x_v,1))
        
        if k<=(size(x_v,2)-1)
            [u_ind]=ind_vort([x_v(l,k),y_v(l,k),z_v(l,k)],[x_v(l,k+1),y_v(l,k+1),z_v(l,k+1)],[x_c(i),y_c(i),z_c(i)],GAMMA);           
            k = k+1;
        else
            l = l+1;
            k = 1;
            [u_ind]=ind_vort([x_v(l,k),y_v(l,k),z_v(l,k)],[x_v(l,k+1),y_v(l,k+1),z_v(l,k+1)],[x_c(i),y_c(i),z_c(i)],GAMMA);  
            k = k+1;
        end 
        
        A(i,j) = dot(u_ind,n);

    end
end

for i = 1:(size(x_c,2)*size(x_c,1))
    b(i,1) = -dot(U_infinity,n);
end
end