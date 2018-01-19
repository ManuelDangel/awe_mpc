% This file sets up some global nmpc parameters
function nmpc = nmpc_init()

% global nmpc  % nmpc is made global for use in simulation

nmpc.N  = 30;        % Horizon Length
nmpc.Ts = 0.1;       % Sample Time

nmpc.x.n = 6;        % State size
nmpc.u.n = 1;        % Input size
nmpc.p.n = 10;       % Parameter size


% Index of states
nmpc.x.index.psi    = 1;  % Azimut Angle
nmpc.x.index.theta  = 2;  % Elevation Angle
nmpc.x.index.gamma  = 3;  % Heading Angle
nmpc.x.index.phi    = 4;  % Roll Angle
nmpc.x.index.vt     = 5;  % Tangential Path Velocity
nmpc.x.index.phi_des= 6;  % Desired (commanded) Roll Angle

% Index of control inputs
nmpc.u.index.dphi   = 1;  % Roll Rate


% Index of parameters
nmpc.p.index.vw                 = 1;  % Wind Velocity
nmpc.p.index.r                  = 2;  % Radius / Tether Length
nmpc.p.index.r_dot              = 3;  % Tether Reel out Speed
nmpc.p.index.circle_azimut      = 4;  % Azimut of Circle Center
nmpc.p.index.circle_elevation   = 5;  % Elevation of Circle Center
nmpc.p.index.circle_angle       = 6;  % Opening Angle of Circle
nmpc.p.index.m                  = 7;  % Mass of the aircraft
nmpc.p.index.clA                = 8;  % Lift coefficient (0.5*rho*A*cl)
nmpc.p.index.cdA                = 9;  % Drag coefficient (0.5*rho*A*cd)
nmpc.p.index.init               =10;  % Initialization boolean in cost function

% Values of parameters
nmpc.p.vw                       = 10;
nmpc.p.r                        = 220;
nmpc.p.r_dot                    = 10*0.23;
nmpc.p.circle_azimut            = 0;
nmpc.p.circle_elevation         = 20/180*pi;
nmpc.p.circle_angle             = atan(85/220);
nmpc.p.m                        = 27.53;
nmpc.p.clA                      = 0.9;
nmpc.p.cdA                      = 0.07;
nmpc.p.init                     = 1;

end