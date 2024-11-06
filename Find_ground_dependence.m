%% INPUT:
clc,clear,close

i = 1;

for h_wing = 0.2:0.2:2
    %---------- WING GEOMETRY: ----------
    c_wing = 3; %chord
    b_wing = 21; %wing span
    % h_wing = 3; %height
    N_wing = 20; %Number panels in chord direction
    M_wing = 20; %Number panels in half wing (wing span direction) EVEN!!!
    Sweep_wing = deg2rad(0); %Sweep angle
    h(i) = h_wing;
    
    %---------- TAIL GEOMETRY: ----------
    c_tail = 1.5; %chord
    b_tail = 6; %wing span
    N_tail = 20; %Number panels in chord direction
    M_tail = 20; %Number panels in half wing (wing span direction) EVEN!!!
    Sweep_tail = deg2rad(0); %Sweep angle
    
    %---------- DISTANCE TAIL-WING: ----------
    x_wt = 5; %x-distance between wing-tail;
    y_wt = 0; %distance LE_wing-LE_tail;
    z_wt = 0.5; %height between wing-tail
    
    
    %---------- FLOW: ----------
    alpha = deg2rad(1); %angle of attack
    beta = deg2rad(0); %angle of sideslip
    U = 1; %freestream intensity
    rho = 1.225; %density
    
    %Create U_infinity as vector
    U_infinity = [U*cos(alpha)*cos(beta), -U*sin(beta), U*sin(alpha)*cos(beta)]; %U_infinity vector
    
    
    %% Build geometry:
    
    %------- WING: -------
    [x_wing,y_wing,z_wing, x_v_wing,y_v_wing,z_v_wing, x_c_wing,y_c_wing,z_c_wing, n_wing,X_c_w,Y_c_w,Z_c_w]=geometry(c_wing,b_wing,N_wing,M_wing,Sweep_wing);
    % IMAGE wing:
    x_wing_im = x_wing;
    y_wing_im = y_wing;
    z_wing_im = z_wing - h_wing;
    x_v_wing_im = x_v_wing;
    y_v_wing_im = y_v_wing;
    z_v_wing_im = z_v_wing - h_wing;
    x_c_wing_im = x_c_wing;
    y_c_wing_im = y_c_wing;
    z_c_wing_im = z_c_wing - h_wing;
    n_wing_im = -n_wing;
    
    %Add the height
    z_wing = z_wing+h_wing;
    z_c_wing = z_c_wing+h_wing;
    z_v_wing = z_v_wing+h_wing;
    Z_c_w = Z_c_w+h_wing;
    
    
    %------- TAIL: -------
    [x_tail,y_tail,z_tail, x_v_tail,y_v_tail,z_v_tail, x_c_tail,y_c_tail,z_c_tail, n_tail,X_c_t,Y_c_t,Z_c_t]=geometry(c_tail,b_tail,N_tail,M_tail,Sweep_tail);
    %Add the height and set position of tail
    x_tail = x_tail+x_wt;
    x_c_tail = x_c_tail+x_wt;
    x_v_tail = x_v_tail+x_wt;
    X_c_t = X_c_t+x_wt;
    
    y_tail = y_tail+y_wt;
    y_c_tail = y_c_tail+y_wt;
    y_v_tail = y_v_tail+y_wt;
    Y_c_t = Y_c_t+y_wt;
    
    z_tail = z_tail+h_wing+z_wt;
    z_c_tail = z_c_tail+h_wing+z_wt;
    z_v_tail = z_v_tail+h_wing+z_wt;
    Z_c_t = Z_c_t+h_wing+z_wt;
    
    
    % IMAGES tail:
    x_tail_im = x_tail;
    y_tail_im = y_tail;
    z_tail_im = z_tail - (2*h_wing+2*z_wt);
    x_v_tail_im = x_v_tail;
    y_v_tail_im = y_v_tail;
    z_v_tail_im = z_v_tail - (2*h_wing+2*z_wt);
    x_c_tail_im = x_c_tail;
    y_c_tail_im = y_c_tail;
    z_c_tail_im = z_c_tail - (2*h_wing+2*z_wt);
    n_tail_im = -n_tail;
    
    
    
%     % Plot geometry 3D
%     figure()
%     surf(x_wing, y_wing, z_wing)
%     hold on
%     plot3(x_c_wing, y_c_wing, z_c_wing,'*b')
%     plot3(x_v_wing, y_v_wing, z_v_wing, '*r')
%     
%     surf(x_tail, y_tail, z_tail)
%     plot3(x_c_tail, y_c_tail, z_c_tail,'*b')
%     plot3(x_v_tail, y_v_tail, z_v_tail, '*r')
%     
%     surf(x_wing_im, y_wing_im, z_wing_im)
%     plot3(x_c_wing_im, y_c_wing_im, z_c_wing_im,'*b')
%     plot3(x_v_wing_im, y_v_wing_im, z_v_wing_im, '*r')
%     
%     surf(x_tail_im, y_tail_im, z_tail_im)
%     plot3(x_c_tail_im, y_c_tail_im, z_c_tail_im,'*b')
%     plot3(x_v_tail_im, y_v_tail_im, z_v_tail_im, '*r')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     axis equal
%     hold off
    
    
    
    %% SOLVE GAMMA:
    
    [A_ww,b_w] = scratc_system(x_c_wing,y_c_wing,z_c_wing,x_v_wing,y_v_wing,z_v_wing,n_wing,U_infinity,1);
    [A_ww_im] = scratc_system(x_c_wing,y_c_wing,z_c_wing,x_v_wing_im,y_v_wing_im,z_v_wing_im,n_wing,U_infinity,1);
    A_ww = A_ww-A_ww_im;
    [A_tw] = scratc_system(x_c_wing,y_c_wing,z_c_wing,x_v_tail,y_v_tail,z_v_tail,n_wing,U_infinity,1);
    [A_tw_im] = scratc_system(x_c_wing,y_c_wing,z_c_wing,x_v_tail_im,y_v_tail_im,z_v_tail_im,n_wing,U_infinity,1);
    A_tw = A_tw+A_tw_im;
    [A_tt,b_t] = scratc_system(x_c_tail,y_c_tail,z_c_tail,x_v_tail,y_v_tail,z_v_tail,n_tail,U_infinity,1);
    [A_tt_im] = scratc_system(x_c_tail,y_c_tail,z_c_tail,x_v_tail_im,y_v_tail_im,z_v_tail_im,n_tail,U_infinity,1);
    A_tt = A_tt+A_tt_im;
    [A_wt] = scratc_system(x_c_tail,y_c_tail,z_c_tail,x_v_wing,y_v_wing,z_v_wing,n_tail,U_infinity,1);
    [A_wt_im] = scratc_system(x_c_tail,y_c_tail,z_c_tail,x_v_wing_im,y_v_wing_im,z_v_wing_im,n_tail,U_infinity,1);
    A_wt = A_wt+A_wt_im;
    
    
    %Build up the system:
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
    
    %Compute all the aerdynamic loads:
    [F_w,Moment_w,C_Lw,C_Dw,C_Mw,Cp_w] = aerodynamic_paramiters(x_wing,y_wing,x_v_wing,y_v_wing,z_v_wing,N_wing, M_wing,GAMMA_w,rho,U_infinity,X_c_w,Y_c_w,Z_c_w);
    [F_t,Moment_t,C_Lt,C_Dt,C_Mt,Cp_t] = aerodynamic_paramiters(x_tail,y_tail,x_v_tail,y_v_tail,z_v_tail,N_tail, M_tail,GAMMA_t,rho,U_infinity,X_c_t,Y_c_t,Z_c_t);
    
    L_total_w = F_w(3);
    D_total_w = F_w(1);
    C_L_total_w(i) = L_total_w/(0.5*rho*c_wing*b_wing*norm(U_infinity)^2);
    C_D_total_w(i) = D_total_w/(0.5*rho*c_wing*b_wing*norm(U_infinity)^2);
    C_M_total_w = - Moment_w(2)/(0.5*rho*c_wing^2*b_wing*norm(U_infinity)^2);
    
    L_total_t = F_t(3);
    D_total_t = F_t(1);
    C_L_total_t(i) = L_total_t/(0.5*rho*c_tail*b_tail*norm(U_infinity)^2);
    C_D_total_t(i) = D_total_t/(0.5*rho*c_tail*b_tail*norm(U_infinity)^2);
    C_M_total_t = - Moment_t(2)/(0.5*rho*c_tail^2*b_tail*norm(U_infinity)^2);
    
    i = i+1;
end

%%
figure()
hold on
for i = 1:length(h)
    plot(h,C_L_total_w,'kd--', 'MarkerFaceColor', 'k')    
    plot(h,C_L_total_t,'kd--', 'MarkerFaceColor', 'b')
end
grid on
title('Lift coefficient VS distance from gorund','FontSize', 15)
xlabel('h','FontSize', 10,'fontweight','bold')
ylabel('C_L','FontSize', 10,'fontweight','bold')
legend('C_L wing','C_L tail','fontsize',10)
% saveas(gcf, 'CL ground dependence','png')

figure()
hold on
for i = 1:length(h)
    plot(h,C_D_total_w,'kd--', 'MarkerFaceColor', 'k')    
    plot(h,C_D_total_t,'kd--', 'MarkerFaceColor', 'b')
end
grid on
title('Drag coefficient VS distance from gorund','FontSize', 15)
xlabel('h','FontSize', 10,'fontweight','bold')
ylabel('C_D','FontSize', 10,'fontweight','bold')
legend('C_D wing','C_D tail','fontsize',10)
% saveas(gcf, 'CD ground dependence','png')