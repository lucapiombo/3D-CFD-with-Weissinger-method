function [u_i]=Biot_Savart(x_v1, x_v2,x_c,gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%APPLY BIOT-SAVART LAW to a vortex:

%INPUTS: coordinates of the extrema of the vortex  x_v1 = [x_v(i,j),y_v(i,j),z_v(i,j)]
%and x_v2 = [x_v(i,j+1),y_v(i,j+1),z_v(i,j+1)], cordinates of the centroid 
%x_c = [x_c,y_c,z_c], and the gamma of the vortex

%OUTPUTS: induced velocity in point x_c due to the vortex (u_i)

%Note: x_v1, x_v2, x_c are NOT the x-component, but are the vectors of the
%left and right extrema of the vortex and of the centorid. 
%Therefore the velocity u_i is also a vector!!!

%NOTE: This is applied only to a segment of the filament, you must recall
%it three times to account for the contribution of the entire horseshoe
%vortex!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r0 = x_v2-x_v1; % Define r0 (vector)

r1 = x_c-x_v1; % Define r1 (vector)
r2 = x_c-x_v2; % Define r2 (vector)

%Apply Biot-Savart:
u_i = gamma/(4*pi)*dot(r0,(r1/norm(r1)-r2/norm(r2)))*(cross(r1,r2)/norm(cross(r1,r2))^2);

end
    
    
