function plot3Dgeometry(geom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the geometry of the wing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

px = geom.panels.x;
py = geom.panels.y;
pz = geom.panels.z;

cx = geom.centroids.x;
cy = geom.centroids.y;
cz = geom.centroids.z;

vx = geom.vortices.x;
vy = geom.vortices.y;
vz = geom.vortices.z;

figure()
surf(px, py, pz)
hold on
plot3(cx, cy, cz,'*b')
plot3(vx, vy, vz, '*r')
xlabel('x','FontSize', 10,'fontweight','bold')
ylabel('y','FontSize', 10,'fontweight','bold')
zlabel('z','FontSize', 10,'fontweight','bold')
axis equal
zlim([-2,3])
hold off
end