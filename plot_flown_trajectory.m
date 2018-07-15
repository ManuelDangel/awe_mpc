% Plot Flown Trajectory at the end
% run after simulation!

r = nmpc.p.r;

% Data to draw reference circle
theta_ref = circle_elevation+circle_angle*sqrt(nmpc.p.r/r)*cos([0:0.01*pi:2*pi]);
%psi_ref = circle_azimut+circle_angle*sqrt(nmpc.p.r/r)*cos([0:0.01*pi:2*pi]);
pre_psi_ref = acos((cos(circle_angle*sqrt(nmpc.p.r/r))-sin(circle_elevation)*sin(theta_ref))./(cos(circle_elevation)*cos(theta_ref)));
psi_ref = [pre_psi_ref(1:100)+circle_azimut, -pre_psi_ref(101:201)+circle_azimut];

x_ref = r*cos(psi_ref).*cos(theta_ref);
y_ref = r*sin(psi_ref).*cos(theta_ref);
z_ref = r*sin(theta_ref);


figure(1)
clf
% set(gcf, 'Position', [screensize(3)*0.5 screensize(4)*0.3 screensize(3)*0.5 screensize(4)*0.7]);
surf(r*x_sphere,r*y_sphere,r*z_sphere,'FaceAlpha',0.5,'FaceColor','interp')
hold on
plot3(x_ref,y_ref,z_ref,'r:');
% plot3(x_pos,y_pos,z_pos,'b');
% plot3(x_pos(1),y_pos(1),z_pos(1),'o');
plot3(pos_history_trac(1,:),pos_history_trac(2,:),pos_history_trac(3,:),'b');
plot3(pos_history_p10(1,:),pos_history_p10(2,:),pos_history_p10(3,:),'r');
% plot3(pos_history_p20(1,:),pos_history_p20(2,:),pos_history_p20(3,:),'g');
axis equal
%title('AWE NMPC Prediction Horizon')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
view(60,30)  % Azimut, Elevation of viewpoint
% legend('sphere at r0','reference at r0','reference tracking','power objective weight 10','power objective weight 20')
legend('sphere at r0','reference at r0','reference tracking','power objective weight 10','Location','Northwest')

orient landscape
print('trajectory_comparison','-dpdf','-bestfit')