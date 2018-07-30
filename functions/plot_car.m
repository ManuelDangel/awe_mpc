%close all;
if 1
    figure(1)
    subplot(6,1,1)
    plot(out.x(:,1));
    ylabel('X')
    subplot(6,1,2)
    plot(out.x(:,2));
    ylabel('Y')
    subplot(6,1,3)
    plot(out.x(:,3));
    ylabel('\psi')
    subplot(6,1,4)
    plot(out.x(:,4));
    ylabel('v_x')
    subplot(6,1,5)
    plot(out.x(:,5));
    ylabel('v_y')
    subplot(6,1,6)
    plot(out.x(:,6));
    ylabel('r')
end
if 1
    figure(3)
    subplot(3,1,1)
    plot(out.x(:,7),'b');
    hold on
    plot(out.u(:,1),'r');
    hold off
    ylabel('\delta')
    subplot(3,1,2)
    plot(out.x(:,8),'b');
    hold on
    plot(out.u(:,2),'r');
    hold off
    ylabel('DC')
    subplot(3,1,3)
    plot(out.u(:,3),'b');
    ylabel('Slack')
end
%%
if 1
    scale=0.5;
    monte_carlo_on = 1;

    f1=figure(2);
    if step<2
        f1.Name='Track';
        f1.Units='normalized';
        f1.OuterPosition = [0 0 1 1];
        xlabel('X [m]')
        ylabel('Y [m]')
        
    end
    if ~exist('h_track_middle')
        h_track_middle = plot(track.middle_line_points(:,1), track.middle_line_points(:,2),'s-','Color', [0.8 0.8 0.8]);
    end
    hold on;
    if ~exist('h_track_right')
        h_track_right = plot(track.right_boundary_cones(:,1), track.right_boundary_cones(:,2),'^-g');
    end
    if ~exist('h_track_left')
        h_track_left = plot(track.left_boundary_cones(:,1),track.left_boundary_cones(:,2),'^-r');
    end
        
    V_max = 15;
    V_min= 5;
    n=100; % Number of points around ellipse
    p=0:pi/n:2*pi; % angles around a circle
    ellipses=[];
    %h6=plot(x_array(:,2,1), x_array(:,2,2),'-ok');
    v_colour = max( 0, min(1, ((out.x(1,4)-V_min)/(V_max-V_min)) ));
    h6=plot(out.x(1:2,1), out.x(1:2,2), 'Color', [v_colour v_colour 1-v_colour ], 'LineWidth', 2);
    horizon_m = zeros(N,2);
    horizon_r = zeros(N,4);
    horizon_l = zeros(N,4);

   
    monte_carlo_data(1, : ,:)= mvnrnd(out.x(1,1:6), confidence*reshape(horizon_state_params.sigma(1,1:6,1:6),6,6), 100);
    for i_N=1:N
         
        %plot ellipses
        [eigvec,eigval] = eig(reshape(horizon_state_params.sigma(i_N,1:2,1:2),2,2)); % Compute eigen-stuff
        xy = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
        x = xy(:,1);
        y = xy(:,2);
        %plot(out.x(i_N,1)+confidence*x, out.x(i_N,2)+confidence*y, 'b');
        %ellipses(i_N, :,:)= [out.x(i_N,1)+confidence*x   out.x(i_N,2)+confidence*y];
        ellipses= [ellipses; out.x(i_N,1)+confidence*x   out.x(i_N,2)+confidence*y];
        
        horizon_m(i_N,:) = [horizon_state_params.P_m(i_N,1), horizon_state_params.P_m(i_N,2)];
        horizon_r(i_N,:) =[horizon_state_params.P_r(i_N,1), horizon_state_params.P_r(i_N,2),scale*input.od(i_N,4),scale*input.od(i_N,5)];
        horizon_l(i_N,:) = [horizon_state_params.P_l(i_N,1),horizon_state_params.P_l(i_N,2),-scale*input.od(i_N,7),-scale*input.od(i_N,8)];
        
        if monte_carlo_on
            next_state = dynamics_6states([reshape(monte_carlo_data(i_N, : ,:),100,6) out.x(i_N,7)*ones(100,1) out.x(i_N,8)*ones(100,1)], car)';
            monte_carlo_data(i_N+1,:,:) = reshape(monte_carlo_data(i_N, : ,:),100,6)+ Ts*next_state(:,1:6);
        end
    end
    if monte_carlo_on ==1
        if exist('h_monte_carlo')
            delete(h_monte_carlo)
        end
       monte_carlo_x=reshape(monte_carlo_data(:,:,1), size(monte_carlo_data,1)*size(monte_carlo_data,2),1);
       monte_carlo_y=reshape(monte_carlo_data(:,:,2), size(monte_carlo_data,1)*size(monte_carlo_data,2),1);

       h_monte_carlo=plot(monte_carlo_x, monte_carlo_y,'r.');

    end    
    if exist('h_x_sol')
       delete(h_x_sol)
    end
    if exist('h_x_guess')
       delete(h_x_guess)
    end
    if exist('h_c_point')
        delete(h_c_point)
    end
    if exist('h_r_arrow')
        delete(h_r_arrow)
    end
    if exist('h_l_arrow')
        delete(h_l_arrow)
    end
     if exist('h_ellipses')
        delete(h_ellipses)
    end
    h_x_sol=plot(out.x(:,1), out.x(:,2),'+-b');
    h_x_guess=plot(input.x(:,1), input.x(:,2),'o-b');
    h_c_point=plot(horizon_m(:,1), horizon_m(:,2),'o', 'Color', [0.8 0.8 0.8]);
    h_r_arrow=quiver(horizon_r(:,1),horizon_r(:,2),horizon_r(:,3),horizon_r(:,4),0,'g');
    h_l_arrow=quiver(horizon_l(:,1),horizon_l(:,2),horizon_l(:,3),horizon_l(:,4),0,'r');
    h_ellipses=plot(ellipses(:,1), ellipses(:,2),'.c');
    axis equal
    %hold off
    

    
end