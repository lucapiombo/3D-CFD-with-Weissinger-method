function [u_ind]=ind_vort(x_v1,x_v2,x_c,gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%APPLY BIOT-SAVART LAW TO THE ENTIRE HORSESHOE VORTEX:

%INPUTS: coordinates of the extrema of horseshoe the horseshoe vortex  
%x_v1 = [x_v(i,j),y_v(i,j),z_v(i,j)] and x_v2 = [x_v(i,j+1),y_v(i,j+1),z_v(i,j+1)], 
%cordinates of the centroid x_c = [x_c,y_c,z_c], and the gamma of the
%horseshoe vortex

%OUTPUTS: induced velocity in point x_c due to the horseshoe vortex (u_i)

%Note: x_v1, x_v2, x_c are NOT the x-component, but are the vectors of the
%left and right extrema of the horseshoe vortex and of the centorid. 
%Therefore the velocity u_i is also a vector!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set extrema of the horshoe vortex at infinite by assuming an high number
x_inf_1 = x_v1+[1*10^6 0 0]; 
x_inf_2 = x_v2+[1*10^6 0 0];

%Recall Biot-Savart law three times:
u_ind = Biot_Savart(x_v1, x_v2,x_c,gamma) + Biot_Savart(x_inf_1, x_v1,x_c,gamma) + Biot_Savart(x_v2,x_inf_2, x_c,gamma);

end