function [u_ind]=ind_vort(v1,v2,centroid,gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%APPLY BIOT-SAVART LAW TO THE ENTIRE HORSESHOE VORTEX:

%INPUTS: coordinates of the extrema of horseshoe the horseshoe vortex  
%v1 = [x_v(i,j),y_v(i,j),z_v(i,j)] and v2 = [x_v(i,j+1),y_v(i,j+1),z_v(i,j+1)], 
%cordinates of the centroid centroid = [centroid,y_c,z_v_inf1c], and the gamma of the
%horseshoe vortex

%OUTPUTS: induced velocity in point centroid due to the horseshoe vortex (u_i)

%Note: v1, v2, centroid are NOT the x-component, but are the vectors of the
%left and right extrema of the horseshoe vortex and of the centorid. 
%Therefore the velocity u_i is also a vector!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set extrema of the horshoe vortex at infinite by assuming an high number
v1_inf = v1+[1*10^6 0 0]; 
v2_inf = v2+[1*10^6 0 0];

%Recall Biot-Savart law three times:
u_ind = biotSavart(v1, v2,centroid,gamma) + ...
        biotSavart(v1_inf, v1,centroid,gamma) + ...
        biotSavart(v2,v2_inf, centroid,gamma);

end