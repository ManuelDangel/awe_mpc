%% AWE NMPC SIMULATION

% Simulation Setup
N=30;
Ts=0.1;

T_Simulation = 10;
stepwise_init = 0;

% Reference Circle
% R=20;
R = atan(90/220);%25/180*pi;% Reference opening angle of the circle

% OnlineData (Parameters)
%vt = 20;
vw = -5;
r = 220;
circle_azimut = -pi/2;
circle_elevation = 20/180*pi;
if stepwise_init
    init = 0;
else
    init = 1;
end


% Data to draw reference circle
psi_ref = circle_azimut+R*cos([0:0.01*pi:2*pi]);
theta_ref = circle_elevation+R*sin([0:0.01*pi:2*pi]);

x_ref = r*cos(psi_ref).*cos(theta_ref);
y_ref = r*sin(psi_ref).*cos(theta_ref);
z_ref = r*sin(theta_ref);


% Initial Positon
% X0 = [psi theta gamma phi vt]
% X0 = [-pi/2+0.3 ,R+circle_elevation+0.1 ,pi/4 ,0.0, 20 ];
X0 = [-pi/2 ,R+circle_elevation ,-pi/2 ,0.0, 20 ];



%% Prepare Simulation
% Fill input vector
input.x0 = X0;
Xref = X0;
input.x = repmat(Xref,N+1,1);%+randn(N+1,4);
%input.x(:,4) = 0.4;
%input.x(:,1) = -pi/2+R*cos(linspace(0,2*pi,N+1)');
%input.x(:,2) = circle_elevation+R*sin(linspace(0,2*pi,N+1)');

input.od = repmat([vw,r,circle_azimut,circle_elevation,init],N+1,1);
%input.od = [repmat([vt,r,circle_azimut,circle_elevation],N/2+1,1);repmat([vt,r,circle_azimut+0.3,circle_elevation],N/2,1)];

%Uref = randn(N,1)*0.1;
Uref = zeros(N,1);
input.u = Uref;

if stepwise_init
    Yref = ones(N,1)*[0,1];
else
    Yref = ones(N,1)*[R,1];
end

input.y = [ Yref, Uref];
input.yN = [R,1];

input.W = diag([1 0 0.01]);
if stepwise_init
    input.WN = diag([0 0]);
else
    input.WN = diag([1 0]);
end


% Constraints
% Control_ub=[max_roll_rate]';
% Control_lb=[-max_roll_rate]';
% 
% State_ub=[ 20  20  10  max_roll_angle*2]';
% State_lb=[-20 -20 -10 -max_roll_angle*2]';
% 
% input.ubValues = repmat(Control_ub,N,1);
% input.lbValues = repmat(Control_lb,N,1);
% 
% input.ubAValues = repmat(State_ub,N,1);
% input.lbAValues = repmat(State_lb,N,1);

[x_sphere,y_sphere,z_sphere] = sphere;
x_sphere = x_sphere(11:21,:);
y_sphere = y_sphere(11:21,:);
z_sphere = z_sphere(11:21,:);



kktValue = [];
objValue = [];

%% Initialization steps
for i=1:N+20
    
    if stepwise_init
        input.od(min(i,N+1),5) = 1;
        input.y(min(i,N),1) = R;
        if i==N
            input.WN = diag([1 0]);
        end
    end

    tic
    output = awe_MPCstep(input); % Solve NMPC
    toc
    
    % Initialize with last output
    input.x = output.x;
    input.u = output.u;
    
    plottingfun(output,x_ref,y_ref,z_ref,x_sphere,y_sphere,z_sphere,r,N);  % plot current trajectory
    
    kktValue = [kktValue output.info.kktValue];
    objValue = [objValue output.info.objValue];
    
    if (i>10 && stepwise_init==0) || i>N+5
    if (abs(objValue(end)-objValue(end-1)) < 10e-2 && kktValue(end) < 10e-5)
        break
    end
    end
    
    pause(Ts-toc)
    
end

%% Run Simulation
cputime=[];
for t=0:Ts:T_Simulation  % Simulation
    
    tic
    output = awe_MPCstep(input); % Solve NMPC
    toc
    
    if false
        input.x = output.x;
        input.u = output.u;
        output = awe_MPCstep(input);
    end

    % Propagate State
    input.x0 = output.x(2,:)+[0,0,0,randn(1,1)*0.0, 0]+randn(1,5)*0.0;
    
    % Propagate Horizon
    input.x(1:end-1,:) = output.x(2:end,:);
    input.x(end,:) = output.x(end,:);
    input.u(1:end-1,:) = output.u(2:end,:);
    input.u(end,:) = output.u(end,:);
    
    plottingfun(output,x_ref,y_ref,z_ref,x_sphere,y_sphere,z_sphere,r,N);  % plot current trajectory
    
    % Logging
    cputime = [cputime,output.info.cpuTime];
    output.info
    %pause(Ts-toc)
end

%% Plot CPU Time
figure(3)
histogram(cputime)

disp('------------------------------------')
disp(['Mean Computation Time: ',num2str(mean(cputime))])
disp(['Max  Computation Time: ',num2str(max(cputime))])
disp(['Min  Computation Time: ',num2str(min(cputime))])
disp('------------------------------------')

% Plotting

% x_out = output.x(:,1);
% y_out = output.x(:,2);
% psi_out = output.x(:,3);
% phi_out = output.x(:,4);
% dphi_out = output.u(:,1);



function plottingfun(output,x_ref,y_ref,z_ref,x_sphere,y_sphere,z_sphere,r,N)
    % plots current trajectory
    
    x_pos = r*cos(output.x(:,1)).*cos(output.x(:,2));
    y_pos = r*sin(output.x(:,1)).*cos(output.x(:,2));
    z_pos = r*sin(output.x(:,2));
    
    screensize=get(0,'Screensize');
    
    figure(1)
    clf
    set(gcf, 'Position', [screensize(3)*0.5 screensize(4)*0.3 screensize(3)*0.5 screensize(4)*0.7]);
    surf(r*x_sphere,r*y_sphere,r*z_sphere,'FaceAlpha',0.5,'FaceColor','interp')
    % plot(cos([0:0.01*pi:2*pi])*R,sin([0:0.01*pi:2*pi])*R,'r')
    hold on
    % plot(output.x(:,1),output.x(:,2),'b')
    % plot(output.x(1,1),output.x(1,2),'o')
    plot3(x_ref,y_ref,z_ref,'r');
    plot3(x_pos,y_pos,z_pos,'b');
    plot3(x_pos(1),y_pos(1),z_pos(1),'o');
    %trajectory_MFILE(x_pos,y_pos,z_pos,zeros(size(x_pos)),zeros(size(x_pos)),zeros(size(x_pos)),1,0)
    %trajectory_MFILE(x_pos(1),y_pos(1),z_pos(1),0,0,0,1,0)
    axis equal
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    % axis([-1.5*R 1.5*R -1.5*R 1.5*R])
    
    figure(2)
    set(gcf, 'Position', [1 screensize(4)*0.3 screensize(3)*0.5 screensize(4)*0.7]);
    subplot(3,1,1); 
    stairs(output.u(:,1)*180/pi);
    grid on; 
    title('Roll Rate [�/s]');
    xlim([0,N+1])
    subplot(3,1,2); 
    stairs(output.x(:,4)*180/pi);
    grid on; 
    title('Roll Angle [�]');
    xlim([0,N+1])
    subplot(3,1,3); 
    stairs(output.x(:,5));
    grid on; 
    title('V tangetial [m/s]');
    xlim([0,N+1])
    
    %view(82.50,2);

end