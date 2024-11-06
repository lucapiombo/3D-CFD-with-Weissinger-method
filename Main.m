clc,clear,close

%% INPUT:
%---------- WING GEOMETRY: ----------
c = 3; %chord
b_span = 18; %wing span
N = 20; %Number panels in chord direction
M = 20; %Number panels in half wing (wing span direction) EVEN!!!
Sweep = deg2rad(0); %Sweep angle


%---------- FLOW: ----------
alpha = deg2rad(1); %angle of attack
beta = deg2rad(0); %angle of sideslip
U = 1; %freestream intensity
rho = 1; %density

%Create U_infinity as vector
U_infinity = [U*cos(alpha)*cos(beta), -U*sin(beta), U*sin(alpha)*cos(beta)]; %U_infinity vector


%% Build geometry:

%------- WING: -------
[x,y,z, x_v,y_v,z_v, x_c,y_c,z_c, n,X_c,Y_c,Z_c]=geometry(c,b_span,N,M,Sweep);


% Plot geometry 3D
figure()
surf(x, y, z)
hold on
plot3(x_c, y_c, z_c,'*b')
plot3(x_v, y_v, z_v, '*r')
xlabel('x','FontSize', 10,'fontweight','bold')
ylabel('y','FontSize', 10,'fontweight','bold')
zlabel('z','FontSize', 10,'fontweight','bold')
axis equal
zlim([-2,3])
hold off




%% SOLVE GAMMA:

[A,b] = scratc_system(x_c,y_c,z_c,x_v,y_v,z_v,n,U_infinity,1);
GAMMA = A\b;
GAMMA = reshape(GAMMA,[2*M,N])';

%% POST PROCESSING:

%Compute all the aerdynamic loads:
[F,Moment,C_L,C_D,C_M,Cp] = aerodynamic_paramiters(x,y,x_v,y_v,z_v,N, M,GAMMA,rho,U_infinity,X_c,Y_c,Z_c);

L_total = F(3);
D_total = F(1);
C_L_total = L_total/(0.5*rho*c*b_span*norm(U_infinity)^2);
C_D_total = D_total/(0.5*rho*c*b_span*norm(U_infinity)^2);
C_M_total = - Moment(2)/(0.5*rho*c^2*b_span*norm(U_infinity)^2);


%FIGURES: 

%---------------- Cp Spanwise direction -----------------
%WING:
figure()
hold on
title('C_p distribution','FontSize', 15)
% for i = 1:N
    plot(y_c(round(N/4),:),Cp(round(N/4),:), '-ob') %if only at 1/4 Chord use only i = round(N/4)
% end
xlabel('Span','FontSize', 10,'fontweight','bold')
ylabel('C_p','FontSize', 10,'fontweight','bold')
grid on
hold off
% saveas(gcf, 'Cp single wing','png')

%---------------- Cp Chordwise direction -----------------
%WING:
figure()
hold on
title('C_p distribution','FontSize', 15)
for i = 1:M
    plot(x_c(:,i),Cp(:,i))
end
grid on
xlabel('Chord','FontSize', 10,'fontweight','bold')
ylabel('C_p','FontSize', 10,'fontweight','bold')
hold off



%---------------- C_L Spanwise direction -----------------
figure()
hold on
% for i = 1:N
    plot(y_c(6,:),C_L(6,:),'-o')
% end
grid on
xlabel('Span','FontSize', 10,'fontweight','bold')
ylabel('C_L','FontSize', 10,'fontweight','bold')
title('C_L distribution along the span','FontSize', 15)
% saveas(gcf, 'CL single wing','png')


%---------------- C_D Spanwise direction -----------------
figure()
hold on
for i = 1:N
    plot(y_c(i,:),C_D(i,:))
end
grid on
xlabel('Span','FontSize', 10,'fontweight','bold')
ylabel('C_D','FontSize', 10,'fontweight','bold')
title('C_D distribution along the span','FontSize', 15)
% saveas(gcf, 'CD single wing','png')



%---------------- GAMMA distribution: -----------------
%COLOR PANELS:
figure()
surf(x,y,z,GAMMA)
hold on
colorbar
title('GAMMA','FontSize', 15)
xlabel('x','FontSize', 10,'fontweight','bold')
ylabel('y','FontSize', 10,'fontweight','bold')
zlabel('z','FontSize', 10,'fontweight','bold')
axis equal
zlim([-2,3])
hold off
% saveas(gcf, 'GAMMA wing','png')


%Show GAMMA in spanwise:
figure()
hold on
for i = 1:N
    plot(y_c(i,:),GAMMA(i,:),'-o')
end
grid on
xlabel('Span','FontSize', 10,'fontweight','bold')
ylabel('GAMMA','FontSize', 10,'fontweight','bold')
title('Gamma distribution along the span','FontSize', 15)


%% VALIDATION (uncommend when you have XFLR5 data):
% Gamma_data = load('GG.txt');
% 
% figure()
% hold on
% plot(y_c(1,:),GAMMA(1,:),'-b')
% plot(Gamma_data(:,1),Gamma_data(:,2),'kd', 'MarkerFaceColor', 'g')
% grid on
% legend('Code','Reference','fontsize',10)
% xlabel('Span','FontSize', 10,'fontweight','bold')
% ylabel('Gamma','FontSize', 10,'fontweight','bold')
% title('Gamma distribution along the span','FontSize', 15)
% hold off
% % saveas(gcf, 'GAMMA wing validation','png')
% 
% 
% 
% Cp_data = load('Cp_wing.txt');
% figure()
% hold on
% plot(y_c(1,:),Cp(1,:),'-k')
% plot(y_c(round(N/4),:),Cp(round(N/4),:),'-b')
% plot(y_c(round(3*N/4),:),Cp(round(3*N/4),:),'-r')
% plot(Cp_data(27:39,1),Cp_data(27:39,2),'kd', 'MarkerFaceColor', 'g')
% plot(Cp_data(1:13,1),Cp_data(1:13,2),'kd', 'MarkerFaceColor', 'g')
% plot(Cp_data(14:26,1),Cp_data(14:26,2),'kd', 'MarkerFaceColor', 'g')
% grid on
% lg = legend('Leading edge','1/4 chord','3/4 chord','Reference','fontsize',10);
% set(lg,'location','best')
% xlabel('Span','FontSize', 10,'fontweight','bold')
% ylabel('C_p','FontSize', 10,'fontweight','bold')
% title('Pressure distribution along the span','FontSize', 15)
% hold off
% %saveas(gcf, 'Cp wing validation','png')
% 
% Cl_data = load('CL.txt');
% figure()
% hold on
% plot(y_c(1,:),C_L(5,:),'-k')
% plot(Cl_data(:,1),Cl_data(:,2),'kd', 'MarkerFaceColor', 'g')
% grid on
% lg = legend('Code','Reference','fontsize',10);
% set(lg,'location','best')
% xlabel('Span','FontSize', 10,'fontweight','bold')
% ylabel('C_L','FontSize', 10,'fontweight','bold')
% title('Lift coefficient distribution along the span','FontSize', 15)
% hold off
% % saveas(gcf, 'Cp wing validation','png')
