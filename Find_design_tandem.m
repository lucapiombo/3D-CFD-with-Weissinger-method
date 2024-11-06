%% INPUT:
clc,clear,close

i = 1;
for b_tail = 5:10
% for S = 0:5:30

    %---------- WING GEOMETRY: ----------
    c_wing = 3; %chord
    b_wing = 21; %wing span
    h_wing = 0; %height
    N_wing = 20; %Number panels in chord direction
    M_wing = 20; %Number panels in half wing (wing span direction) EVEN!!!
    Sweep_wing = deg2rad(20); %Sweep angle
    
    %---------- TAIL GEOMETRY: ----------
    c_tail = 1.5; %chord
%     b_tail = 6; %wing span
    N_tail = 20; %Number panels in chord direction
    M_tail = 20; %Number panels in half wing (wing span direction) EVEN!!!
    Sweep_tail = deg2rad(0); %Sweep angle
%     SWEEP_t(i) = S;
    AR_t(i) = b_tail^2/c_tail;
    
    
    %---------- DISTANCE TAIL-WING: ----------
    x_wt = 5; %x-distance between wing-tail;
    y_wt = 0; %distance LE_wing-LE_tail;
    z_wt = 0.5; %height between wing-tail
%     hz(i) = z_wt;
    
    
    %---------- FLOW: ----------
    alpha = deg2rad(2.6); %angle of attack
    beta = deg2rad(0); %angle of sideslip
    U = 1; %freestream intensity
    rho = 1.225; %density
    
    %Create U_infinity as vector
    U_infinity = [U*cos(alpha)*cos(beta), -U*sin(beta), U*sin(alpha)*cos(beta)]; %U_infinity vector
    
    
    %% Build geometry:
    
    %------- WING: -------
    [x_wing,y_wing,z_wing, x_v_wing,y_v_wing,z_v_wing, x_c_wing,y_c_wing,z_c_wing, n_wing,X_c_w,Y_c_w,Z_c_w]=geometry(c_wing,b_wing,N_wing,M_wing,Sweep_wing);
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
    
    
    
    
    % Plot geometry 3D
%     figure()
%     surf(x_wing, y_wing, z_wing)
%     hold on
%     plot3(x_c_wing, y_c_wing, z_c_wing,'*b')
%     plot3(x_v_wing, y_v_wing, z_v_wing, '*r')
%     surf(x_tail, y_tail, z_tail)
%     plot3(x_c_tail, y_c_tail, z_c_tail,'*b')
%     plot3(x_v_tail, y_v_tail, z_v_tail, '*r')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     axis equal
%     zlim([-4,4])
%     hold off
    
    
    
    %% SOLVE GAMMA:
    
    [A_ww,b_w] = scratc_system(x_c_wing,y_c_wing,z_c_wing,x_v_wing,y_v_wing,z_v_wing,n_wing,U_infinity,1);
    [A_tw] = scratc_system(x_c_wing,y_c_wing,z_c_wing,x_v_tail,y_v_tail,z_v_tail,n_wing,U_infinity,1);
    [A_tt,b_t] = scratc_system(x_c_tail,y_c_tail,z_c_tail,x_v_tail,y_v_tail,z_v_tail,n_tail,U_infinity,1);
    [A_wt] = scratc_system(x_c_tail,y_c_tail,z_c_tail,x_v_wing,y_v_wing,z_v_wing,n_tail,U_infinity,1);
    
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
    
    %Compute all the aerdynamic loads:
    [F_w,Moment_w,C_Lw,C_Dw,C_Mw,Cpw] = aerodynamic_paramiters(x_wing,y_wing,x_v_wing,y_v_wing,z_v_wing,N_wing, M_wing,GAMMA_w,rho,U_infinity,X_c_w,Y_c_w,Z_c_w);
    [F_t,Moment_t,C_Lt,C_Dt,C_Mt,Cpt] = aerodynamic_paramiters(x_tail,y_tail,x_v_tail,y_v_tail,z_v_tail,N_tail, M_tail,GAMMA_t,rho,U_infinity,X_c_t,Y_c_t,Z_c_t);
    
    L_total = F_w(3);
    D_total = F_w(1);
    C_L_total = L_total/(0.5*rho*c_wing*b_wing*norm(U_infinity)^2);
    C_D_total = D_total/(0.5*rho*c_wing*b_wing*norm(U_infinity)^2);
    C_M_total(i) = - Moment_w(2)/(0.5*rho*c_wing^2*b_wing*norm(U_infinity)^2);

    i = i+1;
end


%% POST PROCESSING:



figure()
hold on
for i = 1:length(AR_t)
    plot(AR_t,C_M_total,'ok--', 'MarkerEdgeColor', 'b')
end
grid on
title('Pitching moment VS Aspect Ratio','FontSize', 15)
xlabel('Tail aspect ratio','FontSize', 10,'fontweight','bold')
ylabel('C_{M_{wing}}','FontSize', 10,'fontweight','bold')
% saveas(gcf, 'Moment wing VS AR tail','png')




% figure()
% hold on
% for i = 1:length(SWEEP_t)
%     plot(SWEEP_t,C_M_total,'kd--', 'MarkerFaceColor', 'g')
% end
% grid on
% title('Pitching moment VS Sweep angle','FontSize', 15)
% xlabel('Tail sweep angle','FontSize', 10,'fontweight','bold')
% ylabel('C_{M_{wing}}','FontSize', 10,'fontweight','bold')
% % saveas(gcf, 'Moment wing VS sweep tail','png')







% figure()
% hold on
% plot(y_c_wing(5,:),C_Lw(5,:),'-ok')
% plot(y_c_tail(5,:),C_Lt(5,:),'-ob')
% grid on
% xlabel('Span')
% ylabel('C_L')
% ylim([0,0.12])
% title('C_L distribution along the span')
% 
% 
% % FIGURES:
% 
% %---------------- Cp Spanwise direction -----------------
% 
% %WING:
% figure()
% hold on
% title('C_p Wing')
% for i = 1:N_wing
%     plot(y_c_wing(i,:),Cp_w(i,:)) %if only at 1/4 Chord use only i = round(N/4)
% end
% xlabel('Wing span')
% ylabel('C_p')
% hold off
% 
% %TAIL:
% figure()
% hold on
% title('C_p Tail')
% for i = 1:N_tail
%     plot(y_c_tail(i,:),Cp_t(i,:)) %if only at 1/4 Chord use only i = round(N/4)
% end
% xlabel('Tail span')
% ylabel('C_p')
% hold off
% 
% %WING & TAIL
% figure()
% hold on
% title('C_p Wing VS Tail')
% for i = 1:N_tail
%     plot(y_c_wing(i,:),Cp_w(i,:), 'b') %if only at 1/4 Chord use only i = round(N/4)
%     plot(y_c_tail(i,:),Cp_t(i,:),'r') %if only at 1/4 Chord use only i = round(N/4)
% end
% legend('Wing', 'Tail')
% xlabel('Span')
% ylabel('C_p')
% hold off
% 
% 
% 
% %---------------- Cp Chordwise direction -----------------
% 
% %WING:
% figure()
% hold on
% title('C_p Wing')
% for i = 1:M_wing
%     plot(x_c_wing(:,i),Cp_w(:,i))
% end
% xlabel('Wing chord')
% ylabel('C_p')
% hold off
% 
% %TAIL:
% figure()
% hold on
% title('C_p Tail')
% for i = 11:M_wing
%     plot(x_c_tail(:,i),Cp_t(:,i)) 
% end
% xlabel('Tail chord')
% ylabel('C_p')
% hold off
% 
% %WING & TAIL
% figure()
% hold on
% title('C_p Wing VS Tail')
% for i = 1:1:M_wing
%     plot(x_c_wing(:,i),Cp_w(:,i), 'b')
%     plot(x_c_tail(:,i),Cp_t(:,i),'r')
% end
% legend('Wing', 'Tail')
% xlabel('Chord')
% ylabel('C_p')
% hold off
% 
% 
% %---------------- COLOR PANELS -----------------
% figure()
% surf(x_wing,y_wing,z_wing,GAMMA_w)
% hold on
% surf(x_tail,y_tail,z_tail,GAMMA_t)
% colorbar
% title('GAMMA')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% zlim([-4,4])
% hold off

