%% NMPC SOLVER SETUP

% Setup NMPC Solver using Acado

clear all

acadoSet('problemname', 'awe_mpc');


%% Optimal Control Problem

Ts = 0.1;  % Sampling Time
N = 60;  % Horizon Length
ocp = acado.OCP( 0.0, N*Ts, N );  % Setup Acado problem

DifferentialState x_pos y_pos psi phi
Control dphi

% Parameters
V = 10;
R = 20;

% Differential Equation
f = dot([x_pos; y_pos; psi; phi]) == ...
    [ V*cos(psi); ...
      V*sin(psi); ...
      tan(phi)*9.81/V; ...
      dphi ...
    ];

ocp.setModel(f);

% Constraints
max_roll_angle = 30*pi/180;
max_roll_rate  = max_roll_angle*2;
ocp.subjectTo( -max_roll_angle <= phi  <= max_roll_angle ); % Bounds
ocp.subjectTo( -max_roll_rate  <= dphi <= max_roll_rate  ); % Bounds

% Cost Function
state_output = [sqrt(x_pos*x_pos+y_pos*y_pos)];
control_output = [dphi];

Q_mat = diag([1 1]);
Q = acado.BMatrix(Q_mat);
ocp.minimizeLSQ( Q, [state_output,control_output] );

QN_mat = diag([1]);
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
    cd ..
end
