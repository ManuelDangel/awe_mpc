%% NMPC SOLVER SETUP

% Setup NMPC Solver using Acado

clear all

acadoSet('problemname', 'awe_mpc');


%% Optimal Control Problem

Ts = 0.1;  % Sampling Time
N = 30;  % Horizon Length
ocp = acado.OCP( 0.0, N*Ts, N );  % Setup Acado problem

% DifferentialState x_pos y_pos psi phi
% DifferentialState psi theta gamma phi
DifferentialState psi theta gamma phi vt

Control dphi

% OnlineData vt r circle_azimut circle_elevation init
OnlineData vw r circle_azimut circle_elevation init


% Differential Equation

% 2d plane
% f = dot([x_pos; y_pos; psi; phi]) == ...
%     [ V*cos(psi); ...
%       V*sin(psi); ...
%       tan(phi)*9.81/V; ...
%       dphi ...
%     ];

% Parameters (aircraft)
% rho = 1.225;
% A = 0.39;
% cL = 0.826;
% cD = 0.0862+0.1;
m = 27.53%2.65;
g = 9.81;

clA = 1%0.5*rho*cL*A;
cdA = 0.09%0.5*rho*cD*A;

% % on sphere, fixed radius, fixed velocity
% f = dot([psi; theta; gamma; phi]) == ...
%     [ vt/(r*cos(theta))*sin(gamma); ...
%       vt/r*cos(gamma); ...
%       1/(vt)*(clA*vt*vt/m*sin(phi) + g*cos(theta)*sin(gamma)); ...
%       dphi ...
%     ];

spsi   = sin(psi);
cpsi   = cos(psi);
stheta = sin(theta);
ctheta = cos(theta);
sgamma = sin(gamma);
cgamma = cos(gamma);

% Rotation Matrix from Sphere Local Plane body axis x',y',z' to x,y,z
R_mat = [cpsi -spsi 0 ; spsi cpsi 0 ; 0 0 1]*[-stheta 0 -ctheta ; 0 1 0 ; ctheta 0 -stheta]*[cgamma -sgamma 0 ; sgamma cgamma 0 ; 0 0 1];

%vt_xyz = R_mat*[vt ; 0 ; 0];

% vt_xyz = [ ...  % velocity on sphere in xyz coordinates
%     -r*stheta*cpsi*theta_dot - r*ctheta*spsi*psi_dot; ...
%     -r*stheta*spsi*theta_dot - r*ctheta*cpsi*psi_dot; ...
%     r*ctheta*theta_dot];

vw_xyz = [0;vw;0];  % wind velocity in xyz coordinates

vw_local = R_mat'*vw_xyz;
v_local = [vt ; 0 ; 0] - vw_local;
v = sqrt(v_local'*v_local);
epsilon = asin(v_local(3)/v);
% v_xyz = vt_xyz-vw_xyz;
% v = sqrt(v_xyz'*v_xyz);
% n = [cpsi*ctheta ; spsi*ctheta ; stheta];
% epsilon = asin(v_xyz'*n/v)%v_xyz'*vt_xyz / (v*vt);

% on sphere, fixed radius, velocity is a state
f = dot([psi; theta; gamma; phi; vt]) == ...
    [ vt/(r*ctheta)*sgamma; ...
      vt/r*cgamma; ...
      1/(vt)*(clA*v*v/m*sin(phi) + g*ctheta*sin(gamma)); ...
      dphi; ...
      v*v/m * (clA*sin(epsilon)*cos(phi) - cdA*cos(epsilon) ) - cos(gamma)*ctheta*g% + 5
    ];

ocp.setModel(f);

% Constraints
max_roll_angle = 70*pi/180;
max_roll_rate  = max_roll_angle*2;
max_theta      = 70*pi/180;
min_theta      = -70*pi/180
ocp.subjectTo( -max_roll_angle <= phi   <= max_roll_angle ); % Bounds
ocp.subjectTo( -max_roll_rate  <= dphi  <= max_roll_rate  ); % Bounds
ocp.subjectTo(  min_theta      <= theta <= max_theta      ); % Bounds to prevent gimbal lock
ocp.subjectTo(  10             <= vt    <= 200             ); % Bounds to prevent under/overspeed


% Cost Function
% state_output = [sqrt(x_pos*x_pos+y_pos*y_pos)];
state_output = [...
    sqrt((psi-circle_azimut)*(psi-circle_azimut)+(theta-circle_elevation)*(theta-circle_elevation))*init...
    cos(phi)*cos(theta)*cos(psi-circle_azimut)*clA*v*v...
    ];
control_output = [dphi];

Q_mat = diag([1 0 1]);
Q = acado.BMatrix(Q_mat);
ocp.minimizeLSQ( Q, [state_output,control_output] );

QN_mat = diag([1 0]);
QN = acado.BMatrix(QN_mat);
ocp.minimizeLSQEndTerm( QN, state_output );

% Setup Solver
mpc = acado.OCPexport(ocp)

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


%% Export MPC Solver

if true
    mpc.exportCode( 'export_awe_mpc' );
    copyfile('../../../ACADOtoolkit/external_packages/qpoases', 'export_awe_mpc/qpoases', 'f')
    
    cd export_awe_mpc
    make_acado_solver('../awe_MPCstep')
    make_sfunction('../awe_sfunction')
    cd ..
end
