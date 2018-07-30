%% Logging 
disp('logging')
%% info
info_array_legend= ['status' 'cpuTime' 'kktValue' 'objValue' 'QP_iter' 'QP_violation'];

info_array(step,1) =out.info.status;
info_array(step,2) = out.info.cpuTime;
info_array(step,3) = out.info.kktValue;
info_array(step,4) = out.info.objValue;
info_array(step,5) =out.info.QP_iter;
info_array(step,6) = out.info.QP_violation; 
%% states

input_array(step)=input;

%% Cost and contraints



A_M_l = input.od(:,1);
B_M_l = input.od(:,2);
C_M_l = input.od(:,3);
OD1_l = input.od(:,10);
OD2_l = input.od(:,11);
OD3_l = input.od(:,12);
OD4_l = input.od(:,13);
OD5_l = input.od(:,14);
OD6_l = input.od(:,15);
OD7_l = input.od(:,16);
OD8_l = input.od(:,17);

[derivatives, intermediate]= dynamics_6states(out.x,car);
x_l = out.x(:,1);
y_l = out.x(:,2);
x_dot_l = derivatives(1,:)';
y_dot_l = derivatives(2,:)';


L_contour_l = (  sqrt(sqrt( 0.0000000001+ (OD1_l.*x_l+OD2_l.*x_dot_l*Ts+ OD3_l.*y_l + OD4_l.*y_dot_l*Ts+ OD5_l ).^2+(OD3_l.*x_l+OD4_l.*x_dot_l*Ts + OD6_l.*y_l+ OD7_l.*y_dot_l*Ts+OD8_l ).^2)));
d_middle_l= ((A_M_l.*x_l+B_M_l.*y_l+C_M_l)); %TODO do absoloute value, maybe?
delta_opt_l = input.x-out.x;

%friction_ellipse_r_l  = F_x_l*0.5*1.25 ;  %<=0.000001*(1.1*D_r)^2    
%friction_ellipse_f_l = F_x_l*0.5*1.25 ;  % <= 0.000001*(1.1*D_f)^2
info_array_legend= ['alpha_f' 'alpha_r' 'alpha_c' 'F_fy' 'F_ry' 'F_x' 'L_contour_l' 'd_middle_l'];

data_array(step,:,:)= [intermediate L_contour_l d_middle_l];

delta_opt(step,:,:) = delta_opt_l;
max_L_contour(step)=max(L_contour_l);
