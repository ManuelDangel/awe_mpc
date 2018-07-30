%close all;

%% States
if 1
figure(11)
subplot(6,1,1)
plot(x_array(:,1,1));
ylabel('X')
subplot(6,1,2)
plot(x_array(:,1,2));
ylabel('Y')
subplot(6,1,3)
plot(x_array(:,1,3));
ylabel('\psi')
subplot(6,1,4)
plot(x_array(:,1,4));
ylabel('v_x')
subplot(6,1,5)
plot(x_array(:,1,5));
ylabel('v_y')
subplot(6,1,6)
plot(x_array(:,1,6));
ylabel('r')
end
%%  Inputs
if 1
figure(12)
subplot(3,1,1)
plot(x_array(:,1,7),'b');
hold on
plot(u_array(:,1,1),'r');
hold off
ylabel('\delta')

subplot(3,1,2)
plot(x_array(:,1,8),'b');
hold on
plot(u_array(:,1,2),'r');
hold off
ylabel('DC')

subplot(3,1,3)
plot(sum(u_array(:,:,3),2),'b');
ylabel('Slack')
end
%% Info
if 1
figure(13)
subplot(6,1,1)
plot(info_array(:,1),'b');
ylabel('Status')
subplot(6,1,2)
plot(info_array(:,2),'b');
ylabel('Time')
subplot(6,1,3)
plot(info_array(:,3),'b');
ylabel('KKT')

subplot(6,1,4)
plot(info_array(:,4),'b');
ylabel('Obj')

subplot(6,1,5)
plot(info_array(:,5),'b');
ylabel('itet')

subplot(6,1,6)
plot(info_array(:,6),'b');
ylabel('Viol')
xlabel('Step')
end
%%
if 1
scale=0.5;
figure(14);

%f1.Name='Track';
%f1.Units='normalized';
%f1.OuterPosition = [0 0 1 1];
xlabel('X [m]')
ylabel('Y [m]')

plot(track.middle_line_points(:,1), track.middle_line_points(:,2),'s-', 'Color', [0.8 0.8 0.8])
hold on
plot(track.right_boundary_cones(:,1), track.right_boundary_cones(:,2),'^-g')
plot(track.left_boundary_cones(:,1),track.left_boundary_cones(:,2),'^-r')

plot(x_array(:,1,1),x_array(:,1,2),'-k', 'LineWidth', 2)
hold off
end
%%

figure(15)
boxplot(sqrt(delta_opt(10:end,:,1).^2+delta_opt(10:end,:,2).^2))
grid on
xlabel('Deviation from initial guess [m]')
ylabel('Horizon Points')
