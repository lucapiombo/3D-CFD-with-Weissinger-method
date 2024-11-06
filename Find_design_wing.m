clc,clear,close

i = 1;

% for S = 0:4:32
for b_span = 6:3:30
    %% INPUT:
    
    %---------- WING GEOMETRY: ----------
    c = 3; %chord
%     b_span = 21; %wing span
    N = 20; %Number panels in chord direction
    M = 20; %Number panels in half wing (wing span direction) EVEN!!!
    Sweep = deg2rad(0); %Sweep angle
%     SWEEP(i) = S;
    AR(i) = b_span^2/c;
    
    %---------- FLOW: ----------
    alpha = deg2rad(2.6); %angle of attack
    beta = deg2rad(0); %angle of sideslip
    U = 1; %freestream intensity
    rho = 1.225; %density
    
    %Create U_infinity as vector
    U_infinity = [U*cos(alpha)*cos(beta), -U*sin(beta), U*sin(alpha)*cos(beta)]; %U_infinity vector
    
    
    %% Build geometry:
    
    %------- WING: -------
    [x,y,z, x_v,y_v,z_v, x_c,y_c,z_c, n,X_c,Y_c,Z_c]=geometry(c,b_span,N,M,Sweep);
    
    
    %         % Plot geometry 3D
    %         figure()
    %         surf(x, y, z)
    %         hold on
    %         plot3(x_c, y_c, z_c,'*b')
    %         plot3(x_v, y_v, z_v, '*r')
    %         xlabel('x')
    %         ylabel('y')
    %         zlabel('z')
    %         axis equal
    %         zlim([-2,3])
    %         hold off
    
    
    
    
    %% SOLVE GAMMA:
    
    [A,b] = scratc_system(x_c,y_c,z_c,x_v,y_v,z_v,n,U_infinity,1);
    GAMMA = A\b;
    GAMMA = reshape(GAMMA',[2*M,N])';
    
    %% POST PROCESSING:
    
    %Compute all the aerdynamic loads:
    [F,Moment,C_L,C_D,C_M,Cp] = aerodynamic_paramiters(x,y,x_v,y_v,z_v,N, M,GAMMA,rho,U_infinity,X_c,Y_c,Z_c);
    
    L_total = F(3);
    D_total = F(1);
    C_L_total(i) = L_total/(0.5*rho*c*b_span*norm(U_infinity)^2);
    C_D_total(i) = D_total/(0.5*rho*c*b_span*norm(U_infinity)^2);
    C_M_total = - Moment(2)/(0.5*rho*c^2*b_span*norm(U_infinity)^2);
    
    
    i = i+1;
end
%%
figure()
hold on
for i = 1:length(AR)
    plot(AR,C_D_total,'ok--', 'MarkerEdgeColor', 'b')
end
grid on
title('Induced Drag VS Aspect Ratio','FontSize', 15)
xlabel('Aspect ratio','FontSize', 10,'fontweight','bold')
ylabel('C_{D_i}','FontSize', 10,'fontweight','bold')
% saveas(gcf, 'Induced drag VS AR','png')

% figure()
% hold on
% for i = 1:length(SWEEP)
%     plot(SWEEP,C_D_total,'kd--', 'MarkerFaceColor', 'g')
% end
% grid on
% title('Induced Drag VS Sweep Angle','FontSize', 15)
% xlabel('Sweep angle','FontSize', 10,'fontweight','bold')
% ylabel('C_{D_i}','FontSize', 10,'fontweight','bold')
% saveas(gcf, 'Induced drag VS sweep','png')

% %FIGURES: 
% 
% %---------------- Cp Spanwise direction -----------------
% %WING:
% figure()
% hold on
% title('C_p distribution')
% for i = 1:N
%     plot(y_c(i,:),Cp(i,:)) %if only at 1/4 Chord use only i = round(N/4)
% end
% xlabel('Span')
% ylabel('C_p')
% hold off
% 
% %---------------- Cp Chordwise direction -----------------
% %WING:
% figure()
% hold on
% title('C_p distribution')
% for i = 1:M
%     plot(x_c(:,i),Cp(:,i))
% end
% xlabel('Chord')
% ylabel('C_p')
% hold off
% 
% 
% 
% %---------------- C_L Spanwise direction -----------------
% figure()
% hold on
% for i = 1:N
%     plot(y_c(i,:),C_L(i,:),'-o')
% end
% grid on
% xlabel('Span')
% ylabel('C_L')
% title('C_L distribution along the span')
% 
% 
% %---------------- C_D Spanwise direction -----------------
% figure()
% hold on
% for i = 1:N
%     plot(y_c(i,:),C_D(i,:),'-o')
% end
% grid on
% xlabel('Span')
% ylabel('C_D')
% title('C_D distribution along the span')
% 
% 
% 
% %---------------- GAMMA distribution: -----------------
% %COLOR PANELS:
% figure()
% surf(x,y,z,GAMMA)
% hold on
% colorbar
% title('GAMMA')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% zlim([-2,3])
% hold off
% 
% 
% %Show GAMMA in spanwise:
% figure()
% hold on
% for i = 1:N
%     plot(y_c(i,:),GAMMA(i,:),'-o')
% end
% grid on
% xlabel('Span')
% ylabel('GAMMA')
% title('Gamma distribution along the span')
% 
