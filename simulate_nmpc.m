%% AWE NMPC SIMULATION

% Simulation Setup
clear functions  % Clears compiled functions like .mex
% important to not get stuck in NaN's for consequtive script runs
nmpc = nmpc_init();

N=nmpc.N;
Ts=nmpc.Ts;

T_Simulation = 10;
stepwise_init = 0;

% OnlineData (Parameters)
vw                  = nmpc.p.vw;
r                   = nmpc.p.r;
r_dot               = nmpc.p.r_dot;
circle_azimut       = nmpc.p.circle_azimut;
circle_elevation    = nmpc.p.circle_elevation;
circle_angle        = nmpc.p.circle_angle;

if stepwise_init
    init = 0;
else
    init = 1;
end

% Initial Positon
% X0 = [psi theta gamma phi vt]
% X0 = [-pi/2+0.3 ,circle_angle+circle_elevation+0.1 ,pi/4 ,0.0, 20 ];
X0 = zeros(1,nmpc.x.n);
X0(nmpc.x.index.psi)    = circle_azimut;
X0(nmpc.x.index.theta)  = circle_angle+circle_elevation;
X0(nmpc.x.index.gamma)  = -pi/2;
X0(nmpc.x.index.phi)    = 0;
X0(nmpc.x.index.vt)     = 20;




%% Prepare Simulation
% Fill input vector
input.x0 = X0;
Xref = X0;
input.x = repmat(Xref,N+1,1);

%input.od = repmat([vw,r,circle_azimut,circle_elevation,init],N+1,1);
input.od = zeros(N+1,nmpc.p.n)
input.od(:,nmpc.p.index.vw)                 = repmat(vw,N+1,1);
for i=1:N+1
    input.od(i,nmpc.p.index.r)              = r + (i-1)*Ts*r_dot;
end
input.od(:,nmpc.p.index.r_dot)              = repmat(r_dot,N+1,1);
input.od(:,nmpc.p.index.circle_azimut)      = repmat(circle_azimut,N+1,1);
input.od(:,nmpc.p.index.circle_elevation)   = repmat(circle_elevation,N+1,1);
input.od(:,nmpc.p.index.circle_angle)       = circle_angle*sqrt(nmpc.p.r./input.od(:,nmpc.p.index.r));  % repmat(circle_angle,N+1,1);
input.od(:,nmpc.p.index.m)                  = repmat(nmpc.p.m,N+1,1);
input.od(:,nmpc.p.index.clA)                = repmat(nmpc.p.clA,N+1,1);
input.od(:,nmpc.p.index.cdA)                = repmat(nmpc.p.cdA,N+1,1);
input.od(:,nmpc.p.index.init)               = repmat(init,N+1,1);


Uref = zeros(N,nmpc.u.n);
input.u = Uref; 

Yref = ones(N,1)*[0,1];

input.y = [Yref, Uref];
input.yN = [0,1];

input.W = diag([1 0 0.01]);
if stepwise_init
    input.WN = diag([0 0]);
else
    input.WN = diag([1 0]);
end


% % Change last state costs to track trajectory well:
% input.od(end-6:end-1,5)=10;


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
x_sphere = x_sphere(11:21,6:16);
y_sphere = y_sphere(11:21,6:16);
z_sphere = z_sphere(11:21,6:16);



kktValue = [];
objValue = [];

%% Initialization steps
for i=1:N+20
    
    if stepwise_init
        input.od(min(i,N+1),5) = 1;
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
    
    plottingfun(output,x_sphere,y_sphere,z_sphere,circle_azimut,circle_elevation,circle_angle,r,r_dot,nmpc);  % plot current trajectory
    
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

% % Cost Weighting for Power optimization
%input.W = diag([1 0.15 0.01]);

cputime=[];
for t=0:Ts:T_Simulation  % Simulation
    
    tic
    output = awe_MPCstep(input); % Solve NMPC
    toc
    
    % Propagate State
    input.x0 = output.x(2,:);%+[0,0,0,randn(1,1)*0.0, 0]+randn(1,5)*0.0;
    
    % Propagate Horizon
    input.x(1:end-1,:) = output.x(2:end,:);
    input.x(end,:) = output.x(end,:);
    input.u(1:end-1,:) = output.u(2:end,:);
    input.u(end,:) = output.u(end,:);
    
    % Wrap Heading Angle Horizon 
if input.x0(nmpc.x.index.gamma)-input.x(1,nmpc.x.index.gamma) < -pi
    input.x(:,nmpc.x.index.gamma) = input.x(:,nmpc.x.index.gamma) - 2*pi*ones(size(input.x(:,nmpc.x.index.gamma)))
end
if input.x0(nmpc.x.index.gamma)-input.x(1,nmpc.x.index.gamma) > pi
    input.x(:,nmpc.x.index.gamma) = input.x(:,nmpc.x.index.gamma) + 2*pi*ones(size(input.x(:,nmpc.x.index.gamma)))
end
    plottingfun(output,x_sphere,y_sphere,z_sphere,circle_azimut,circle_elevation,circle_angle,r,r_dot,nmpc);  % plot current trajectory
    
    % Propagate Tether Length
    r = r + Ts * r_dot;
    for i=1:N+1
        input.od(i,nmpc.p.index.r)              = r + (i-1)*Ts*r_dot;
    end
    input.od(:,nmpc.p.index.circle_angle)       = circle_angle*sqrt(nmpc.p.r./input.od(:,nmpc.p.index.r));  % repmat(circle_angle,N+1,1);

    
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


function plottingfun(output,x_sphere,y_sphere,z_sphere,circle_azimut,circle_elevation,circle_angle,r,r_dot,nmpc)
    % plots current trajectory
    
    % Data to draw reference circle
    theta_ref = circle_elevation+circle_angle*sqrt(nmpc.p.r/r)*sin([0:0.01*pi:2*pi]);
    psi_ref = circle_azimut+circle_angle*sqrt(nmpc.p.r/r)*cos([0:0.01*pi:2*pi]);

    x_ref = r*cos(psi_ref).*cos(theta_ref);
    y_ref = r*sin(psi_ref).*cos(theta_ref);
    z_ref = r*sin(theta_ref);
    
    r_traj = zeros(nmpc.N+1,1);
    for i=1:nmpc.N+1  % Calculate the expanding trajectory r
        r_traj(i) = r + (i-1)*nmpc.Ts*r_dot;
    end

    x_pos = r_traj.*cos(output.x(:,nmpc.x.index.psi)).*cos(output.x(:,nmpc.x.index.theta));
    y_pos = r_traj.*sin(output.x(:,nmpc.x.index.psi)).*cos(output.x(:,nmpc.x.index.theta));
    z_pos = r_traj.*sin(output.x(:,nmpc.x.index.theta));
    
    screensize=get(0,'Screensize');
    
    figure(1)
    clf
    % set(gcf, 'Position', [screensize(3)*0.5 screensize(4)*0.3 screensize(3)*0.5 screensize(4)*0.7]);
    surf(r*x_sphere,r*y_sphere,r*z_sphere,'FaceAlpha',0.5,'FaceColor','interp')
    % plot(cos([0:0.01*pi:2*pi])*circle_angle,sin([0:0.01*pi:2*pi])*circle_angle,'r')
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
    view(60,30)  % Azimut, Elevation of viewpoint
    
    figure(2)
    % set(gcf, 'Position', [1 screensize(4)*0.3 screensize(3)*0.5 screensize(4)*0.7]);
    subplot(3,1,1); 
    stairs(output.u(:,nmpc.u.index.dphi)*180/pi);
    grid on; 
    title('Roll Rate [�/s]');
    xlim([0,nmpc.N+1])
    subplot(3,1,2); 
    stairs(output.x(:,nmpc.x.index.phi)*180/pi);
    grid on; 
    title('Roll Angle [�]');
    xlim([0,nmpc.N+1])
    subplot(3,1,3); 
    stairs(output.x(:,nmpc.x.index.vt));
    grid on; 
    title('V tangetial [m/s]');
    xlim([0,nmpc.N+1])

end
