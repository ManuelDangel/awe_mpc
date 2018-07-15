%% NMPC SOLVER SETUP

% Setup NMPC Solver using Acado

clear all

acadoSet('problemname', 'awe_mpc');


%% Optimal Control Problem

nmpc = nmpc_init();

Ts = nmpc.Ts;  % Sampling Time
N  = nmpc.N;   % Horizon Length
ocp = acado.OCP( 0.0, N*Ts, N );  % Setup Acado problem

% DifferentialState x_pos y_pos psi phi
% DifferentialState psi theta gamma phi
DifferentialState psi theta gamma phi vt phi_des
  
Control dphi phi_slack theta_slack

% OnlineData vt r circle_azimut circle_elevation init
OnlineData vw r r_dot circle_azimut circle_elevation circle_angle m clA cdA phi_freq wind_azimut thrust_power weight_tracking weight_power


% Parameters (aircraft)
g = 9.81;

% Intermediate States
spsi   = sin(psi);
cpsi   = cos(psi);
stheta = sin(theta);
ctheta = cos(theta);
sgamma = sin(gamma);
cgamma = cos(gamma);

% Rotation Matrix from Sphere Local Plane body axis x',y',z' to x,y,z
R_mat = [cpsi -spsi 0 ; spsi cpsi 0 ; 0 0 1]*[-stheta 0 -ctheta ; 0 1 0 ; ctheta 0 -stheta]*[cgamma -sgamma 0 ; sgamma cgamma 0 ; 0 0 1];

vw_xyz = [vw*cos(wind_azimut);vw*sin(wind_azimut);0];  % wind velocity in xyz coordinates

vw_local = R_mat'*vw_xyz;
v_local = [vt ; 0 ; -r_dot] - vw_local;
v = sqrt(v_local'*v_local);
epsilon = asin(v_local(3)/v);
beta = atan(v_local(2)/vt);

% Differential Equation
% on sphere, fixed radius, 3D wind velocity, path side slip
if nmpc.flip
    f = dot([psi; theta; gamma; phi; vt ; phi_des]) == ...
    [ vt/(r*ctheta)*sgamma; ...
      vt/r*cgamma; ...
      1/(vt)*(clA*v*v/m*sin(phi)*cos(beta) + g*(-cpsi*cgamma - spsi*stheta*sgamma)); ...
      (phi_des*1.0-phi)*phi_freq; ...
      v*v/m * (clA*(sin(epsilon)*cos(phi)*cos(beta) + sin(phi)*sin(beta)) - cdA*cos(epsilon) ) + g*(-cpsi*sgamma + spsi*stheta*cgamma) + thrust_power/vt; ...% + 5 
      dphi
    ];
else
    f = dot([psi; theta; gamma; phi; vt ; phi_des]) == ...
    [ vt/(r*ctheta)*sgamma; ...
      vt/r*cgamma; ...
      1/(vt)*(clA*v*v/m*sin(phi)*cos(beta) + g*ctheta*sgamma); ...
      (phi_des*1.0-phi)*phi_freq; ...
      v*v/m * (clA*(sin(epsilon)*cos(phi)*cos(beta) + sin(phi)*sin(beta)) - cdA*cos(epsilon) ) - cgamma*ctheta*g; ...% + 5 
      dphi
    ];
end


ocp.setModel(f);

% Constraints
max_roll_angle = 40*pi/180;  % Soft Constraint
max_roll_rate  = max_roll_angle*2;
max_theta      = 70*pi/180;
min_theta      = 10*pi/180;  % Has to be chosen well above zero because of quadratic slack soft constraint
if nmpc.flip
   min_theta      = -70*pi/180;
end
ocp.subjectTo( -max_roll_rate  <= dphi  <= max_roll_rate  ); % Bounds on input (hard constraint)
ocp.subjectTo( max_roll_angle  >= phi-phi_slack ); % Upper Bound on roll angle (soft constraint)
ocp.subjectTo( -max_roll_angle <= phi+phi_slack ); % Lower Bound on roll angle (soft constraint)
ocp.subjectTo(  max_theta      >= theta-theta_slack); % Upper Bound on theta to prevent going close to gimbal lock (soft constraint)
ocp.subjectTo(  min_theta      <= theta+theta_slack); % Lower Bound on theta to prevent hitting the ground (soft constraint)
ocp.subjectTo(  -85*pi/180     <= theta <= 85*pi/180      ); % Bounds to prevent gimbal lock
ocp.subjectTo(  5              <= vt                      ); % Bounds to prevent underspeed and divisions by zero
ocp.subjectTo(  0              <= phi_slack               ); % Bounds on Slack
ocp.subjectTo(  0              <= theta_slack             ); % Bounds on Slack


% Cost Function
thrust = thrust_power/vt + 1;  % prevent blowup in power objective if vw is zero.
power_optimum    = clA * (clA/cdA*((vw+thrust)*2/3))^2*(vw+thrust)/3;
power_actual     = (clA*cos(phi)*cos(epsilon) + cdA*sin(epsilon))*v^2*r_dot;
power_potential  = m*g*vt*cgamma*ctheta;
energy_optimum   = m*g*r + 0.5*m*clA/cdA*((vw+thrust))^2; % MISTAKE!!!
energy_potential = m*g*r*stheta + 0.5*m*vt^2;
xyz_center   = [cos(circle_elevation)*cos(circle_azimut), cos(circle_elevation)*sin(circle_azimut), sin(circle_elevation)];
xyz_position = [ctheta*cpsi, ctheta*spsi, stheta];
angle_from_center = acos(xyz_center*xyz_position');
state_output = [...
    (angle_from_center-circle_angle) * weight_tracking...
    sqrt((10*power_optimum - ( power_actual  + 0.0 * power_potential ))/power_optimum ) * weight_power...
    ];
control_output = [dphi, phi_slack, theta_slack];
end_output = [sqrt( (10*energy_optimum - 1.0*energy_potential)/power_optimum) * weight_power ];

Q_mat = diag([1 0 1 1 1]);
Q = acado.BMatrix(Q_mat);
ocp.minimizeLSQ( Q, [state_output,control_output] );

QN_mat = diag([1 0 0]);
QN = acado.BMatrix(QN_mat);
ocp.minimizeLSQEndTerm( QN, [state_output end_output] );

% Setup Solver
mpc = acado.OCPexport(ocp);

mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'INTEGRATOR_TYPE',             'INT_IRK_GL4'       );  % Implicit 4th order RK
mpc.set( 'DYNAMIC_SENSITIVITY',         'FORWARD'           );  % Default and only one that produces reasonable results
mpc.set( 'LINEAR_ALGEBRA_SOLVER',       'GAUSS_LU'          );  % Default and only one that produces reasonable results
mpc.set( 'NUM_INTEGRATOR_STEPS',        N                   );
% mpc.set( 'IMPLICIT_INTEGRATOR_NUM_ITS', 10                  );
mpc.set( 'QP_SOLVER',                   'QP_QPOASES'    	);
mpc.set( 'HOTSTART_QP',                 'NO'             	);

% mpc.set( 'LEVENBERG_MARQUARDT', 		 1e-10				);
% mpc.set( 'CG_HARDCODE_CONSTRAINT_VALUES','NO' 			);  % Default is yes
mpc.set( 'SPARSE_QP_SOLUTION',          'FULL_CONDENSING'   );
% mpc.set( 'CG_USE_OPENMP',             'YES' );  % Paralellization
% mpc.set( 'GENERATE_SIMULINK_INTERFACE', 'YES' );

%% Export MPC Solver

if true
    mpc.exportCode( 'export_awe_mpc' );
    copyfile('../../../ACADOtoolkit/external_packages/qpoases', 'export_awe_mpc/qpoases', 'f')
    
    cd export_awe_mpc
    make_acado_solver('../awe_MPCstep')
    % make_sfunction('../awe_sfunction')
    cd ..
end
