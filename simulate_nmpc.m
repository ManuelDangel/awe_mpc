%% AWE NMPC SIMULATION

% Simulation Setup
T_Simulation = 30;

% Reference Circle
R=20;

% Initial Positon
X0 = [15 ,0 ,pi/2 ,0 ];




%% Prepare Simulation
% Fill input vector
input.x0 = X0;
Xref = X0;
input.x = repmat(Xref,N+1,1);%+randn(N+1,4);
input.od = [];

%Uref = randn(N,1)*0.1;
Uref = zeros(N,1);
input.u = Uref;

Yref = ones(N,1)*R;
input.y = [ Yref, Uref];
input.yN = R;

input.W = diag([1 1]);
input.WN = diag([1]);

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


%% Initialization steps
for i=1:10
    
    tic
    output = awe_MPCstep(input); % Solve NMPC
    toc
    
    % Initialize with last output
    input.x = output.x;
    input.u = output.u;
    
    % Plotting
    figure(1)
    clf
    plot(cos([0:0.01*pi:2*pi])*R,sin([0:0.01*pi:2*pi])*R,'r')
    hold on
    plot(output.x(:,1),output.x(:,2),'b')
    plot(output.x(1,1),output.x(1,2),'o')
    xlabel('x [m]')
    ylabel('y [m]')
    axis([-1.5*R 1.5*R -1.5*R 1.5*R])
    pause(Ts-toc)
    
end

%% Run Simulation
cputime=[];
for t=0:Ts:T_Simulation  % Simulation
    
    tic
    output = awe_MPCstep(input); % Solve NMPC
    toc
    
    % Propagate State
    input.x0 = output.x(2,:)+[0,0,0,randn(1,1)*0]+randn(1,4)*0.0;
    
    % Propagate Horizon
    input.x(1:end-1,:) = output.x(2:end,:);
    input.x(end,:) = output.x(end,:);
    input.u(1:end-1,:) = output.u(2:end,:);
    input.u(end,:) = output.u(end,:);
    
    % Plotting
    figure(1)
    clf
    plot(cos([0:0.01*pi:2*pi])*R,sin([0:0.01*pi:2*pi])*R,'r')
    hold on
    plot(output.x(:,1),output.x(:,2),'b')
    plot(output.x(1,1),output.x(1,2),'o')
    xlabel('x [m]')
    ylabel('y [m]')
    axis([-1.5*R 1.5*R -1.5*R 1.5*R])
    
    % Logging
    cputime = [cputime,output.info.cpuTime];
    
    pause(Ts-toc)
end

%% Plot CPU Time
figure(2)
histogram(cputime)

disp('------------------------------------')
disp(['Mean Computation Time: ',num2str(mean(cputime))])
disp(['Max Computation Time: ',num2str(max(cputime))])
disp(['Min Computation Time: ',num2str(min(cputime))])
disp('------------------------------------')

% Plotting

% x_out = output.x(:,1);
% y_out = output.x(:,2);
% psi_out = output.x(:,3);
% phi_out = output.x(:,4);
% dphi_out = output.u(:,1);
