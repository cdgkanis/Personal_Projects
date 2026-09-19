clc;
clear;
close all;

%% Initialization

% Robot Parameters
Robot_params = struct();
Robot_params.DIM =  2;
Robot_params.m_i = [2 2]; % [-] No. Frequencies for every link i
Robot_params.L   = [1 1]; % [m]
Robot_params.r   = [0.1 0.1]; % [m] Radii of Equivalent Cylinders

% h := Mass and kinematik Inertia of the hubs
% p := Mass and kinematik Inertia of the payload
Robot_params.Mh = [1 1]; %[Kg]
Robot_params.Mp =  5; %[Kg]
Robot_params.Jh = [0.1 0.1]; %[Kgm^2]
Robot_params.Jp =  0.0005; %[Kgm^2]

Robot_params.End_effector.Kx = 1e7;
Robot_params.End_effector.Ky = 1e7;
Robot_params.End_effector.Cx = 1e4;
Robot_params.End_effector.Cy = 1e4;
Robot_params.End_effector.M_pe = 3.0;

Robot_params.Controler.omega = 80;
Robot_params.Controler.zeta = 0.707;

% Material Parameters
Material_params = struct();
Material_params.E = 200e9; %[Pa]
Material_params.ro = 7800; %[Kg/m^3]

%% Simulation setup

mov_params = struct();
mov_params.cx = 0.2;% Path starting position @ X
mov_params.cy = 0.0;% Path starting position @ Y
mov_params.a = 1.0; % Path scaling factor
mov_params.T = 2; % Path traversing time (currently same as the sim duration)
mov_params.bidirectional = false;

run_options = struct();
run_options.live_plots = true;
run_options.animate = false;
run_options.profile = true;
run_options.spectogram = 'relative'; % 'relative' || 'absolute'

%% Construct and initialize the robot model

Robot = Robot(Robot_params, Material_params);

if Robot.N ~= 2
    error('LightWeightRoboticArmSims: the Cartesian reference is defined for a two-link robot.');
end

Robot.Print_properties();
Robot.Calculate_End_Characteristics();
Robot.Print_End_properties();

% Distal links must be solved before proximal links because their mass and
% inertia enter the proximal-link boundary conditions.
for i = Robot.N:-1:1
    Robot.Calculate_Eigen_Frequencies(i);
end
Robot.Print_frequency_properties();

for i = Robot.N:-1:1
    Robot.Find_C_Coeficients(i);
end
Robot.Check();

Robot.Calculate_K_matrix();
Robot.Calculate_D_matrix();

%% Desired Cartesian trajectory and rigid-joint reference

dt = 1/2^13;

t_span = [0 mov_params.T];
N_time = t_span(2)/dt;
t_grid = linspace(t_span(1), t_span(2), ceil(N_time));

[theta_ref, omega_ref, alpha_ref] = simulate_robot_path(Robot, t_grid, 1, mov_params);
t_ref = linspace(t_span(1), t_span(2), size(theta_ref, 2));

interp_q     = griddedInterpolant(t_ref, theta_ref', 'linear', 'nearest');
interp_omega = griddedInterpolant(t_ref, omega_ref', 'linear', 'nearest');
interp_alpha = griddedInterpolant(t_ref, alpha_ref', 'linear', 'nearest');
%% Initial conditions: x = [q_dot; V_dot_E; q; V_E]

Nq = Robot.NDoFs;

q0 = zeros(Nq,1);
q0(1:Robot.N) = theta_ref(:,1);

q_dot0 = zeros(Nq,1);
Robot.setState(q0,q_dot0);
[r3,~] = Robot.FW_Kinematics_2D(q0,q_dot0,2,Robot.L(2));

vx_E0     = r3(1);
vy_E0     = r3(2);
vx_dot_E0 = 0;
vy_dot_E0 = 0;


x0 = [q_dot0;
      vx_dot_E0;
      vy_dot_E0;
      q0;
      vx_E0;
      vy_E0];

Robot.setState(q0, q_dot0);
Robot.Calculate_H_matrix();

Robot.Print_K_matrix();
Robot.Print_D_matrix();
Robot.Print_B_matrix();
Robot.Print_H_matrix();

%% ODE options and optional live plots

if run_options.live_plots
    live = init_live_plots(Robot, theta_ref, t_grid);
    output_function = @(t,x,flag) live_plot_callback( t, x, flag, Robot, live);
else
    output_function = [];
end

options_ODE = odeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-8, ...
    'MaxStep', 1e-3, ...
    'Vectorized', 'off', ...
    'Stats', 'on', ...
    'OutputFcn', output_function);

%% Integrate

if run_options.profile
    profile clear;
    profile on;
end

tSolve = tic;
[Time_out, X] = ode15s(@(t,x) robot_ode( ...
    t, x, Robot, interp_q, interp_omega, interp_alpha), ...
    t_grid, x0, options_ODE);
solveTime = toc(tSolve);

fprintf('Total ODE solve time: %.6f s\n', solveTime);

if run_options.profile
    profile off;
    profile viewer;
end

%% Post-processing

[tip_position, tip_velocity, EEF_position, EEF_velocity, rel_pos, rel_vel] = postprocess_tip_kinematics(Time_out, X, Robot);

plot_cartesian_path(Robot, tip_position, theta_ref);

plot_tip_quantities( Time_out, ...
    tip_position, tip_velocity, ...
    EEF_position, EEF_velocity, ...
    rel_pos, rel_vel);

if run_options.spectogram == "relative"
    plot_relative_stft(Time_out, rel_pos, Robot.f_bar, Robot.f_sd);
    
elseif run_options.spectogram == "absolute"
    plot_absolute_stft(Time_out, EEF_position);
end


[Tau_hist, P_hist] = postprocess_torques( Time_out, X, Robot, interp_q, interp_omega, interp_alpha);
plot_torque_power(Time_out, Tau_hist, P_hist);

if run_options.animate
    animate_robot(Robot, Time_out, X, theta_ref);
end


%% ========================================================================
%  Dynamics and controller
% =========================================================================

function xdot = robot_ode(t, x, Robot, interp_q, interp_omega, interp_alpha)
    
    Nq = Robot.NDoFs;
    Nr = Robot.N;
    
    % State:
    % x = [q_dot; vX_dot_E; vY_dot_E; q; vX_E; vY_E]^T
    q_dot   = x(1:Nq);
    v_dot(:,1) = x(Nq+1:Nq+2);
    q   = x(Nq+3:2*Nq+2);
    v(:,1) = x(2*Nq+3:2*Nq+4);
    
    Robot.setState(q, q_dot);
    Robot.Calculate_H_matrix();
    [B, D, K, H] = Robot.getDynamicMatrices();
    
    
    %Adding the end effector 
    EEF = Robot.End_effector;
    
    [r_EEF,r_EEF_dot,~] = Robot.FW_Kinematics_2D(q, q_dot, 2, Robot.L(2));
    
    K_EEF = [EEF.Kx 0; 0 EEF.Ky];
    C_EEF = [EEF.Cx 0; 0 EEF.Cy];
    
    J = Robot.shape_functions();
    
    K_BASE = -J'*K_EEF;
    C_BASE = -J'*C_EEF;

    r_res    = r_EEF     - J*q;
    rdot_res = r_EEF_dot - J*q_dot;
    
    F_EEF  = -K_EEF*r_res - C_EEF*rdot_res;
    F_BASE =  J'*(K_EEF*r_res + C_EEF*rdot_res);
    
    K_ = [K+J'*K_EEF*J K_BASE;K_BASE' K_EEF];
    C_ = [D+J'*C_EEF*J C_BASE;C_BASE' C_EEF];
    M_ = blkdiag(B, EEF.M_pe * eye(2));
    H_ = [H;0;0] + [F_BASE;F_EEF];

    % eigen frequencies
    lambda = polyeig(K_, C_, M_);
    lambda_pos = lambda(imag(lambda) >= 0);
    omega_d = imag(lambda_pos);
    Robot.renew_statistics(omega_d);
    
    % Full augmented state used by the inverse-dynamics controller.
    % The LQR feedback still acts only on the two rigid joints, but the
    % feedforward term must retain the flexible and EEF coupling terms.
    z     = [q; v];
    z_dot = [q_dot; v_dot];

    Control_input = zeros(Nq+2,1);
    % Controller and current augmented matrices
    tau_joint = compute_control_torque( ...
        t, z, z_dot, Robot, interp_q, interp_omega, interp_alpha, ...
        M_, C_, K_, H_);

    % Control Input
    Control_input(1:Nr) = tau_joint(:);
    Q = Control_input;
    
    z_ddot = M_ \ (Q - C_*z_dot - K_*z - H_);
    
    % State derivative:
    % x = [q_dot; v_dot; q; v]
    % xdot = [q_ddot; v_ddot; q_dot; v_dot]
    xdot = [z_ddot; q_dot; v_dot];
end

function tau_joint = compute_control_torque( ...
        t, z, z_dot, Robot, interp_q, interp_omega, interp_alpha, ...
        M, D, K, H)
    
    q_d = interp_q(t).';
    q_dot_d = interp_omega(t).';
    q_ddot_d = interp_alpha(t).';
    
    omega_n = Robot.Controler.omega;
    zeta = Robot.Controler.zeta;
    Kp = diag(repmat(omega_n^2,1,Robot.N));
    Kd = diag(repmat(2*zeta*omega_n,1,Robot.N));
    
    Q = blkdiag(Kp, Kd);   % 4-by-4 rigid-state weighting
    R = eye(2);            % 2-by-2 control weighting

    % IMPORTANT: pass the complete augmented matrices and state to LQR.
    % LQR uses only the first two coordinates for feedback, while its
    % inverse-dynamics term retains all unactuated flexible/EEF coupling.
    [Tau1, Tau2] = LQR( ...
        M, D, K, H, z, z_dot, ...
        q_d, q_dot_d, q_ddot_d, Q, R);
    tau_joint = [Tau1; Tau2];

end

%% ========================================================================
%  Post-processing and visualization
% =========================================================================

function [tip_position, tip_velocity, ...
          EEF_position, EEF_velocity, ...
          rel_pos, rel_vel] = ...
          postprocess_tip_kinematics(t, X, Robot)

    nt = numel(t);
    Nq = Robot.NDoFs;

    tip_position = zeros(nt,2);
    tip_velocity = zeros(nt,2);

    EEF_position = zeros(nt,2);
    EEF_velocity = zeros(nt,2);

    rel_pos = zeros(nt,2);
    rel_vel = zeros(nt,2);

    for k = 1:nt

        q_dot = X(k,1:Nq).';
        q     = X(k,Nq+3:2*Nq+2).';

        [r3,r3_dot] = Robot.FW_Kinematics_2D( ...
            q, q_dot, Robot.N, Robot.L(Robot.N));

        % Robot tip
        tip_position(k,:) = r3.';
        tip_velocity(k,:) = r3_dot.';

        % Separate EEF mass
        EEF_position(k,:) = X(k,2*Nq+3:2*Nq+4);
        EEF_velocity(k,:) = X(k,Nq+1:Nq+2);

        % Relative spring/damper deformation
        rel_pos(k,:) = EEF_position(k,:) ...
                     - tip_position(k,:);

        rel_vel(k,:) = EEF_velocity(k,:) ...
                     - tip_velocity(k,:);
    end
end

function [tau_history, power_history] = postprocess_torques( t, X, Robot, interp_q, interp_omega, interp_alpha)
    
    nt = numel(t);
    Nq = Robot.NDoFs;
    tau_history = zeros(nt,Robot.N);
    power_history = zeros(nt,Robot.N);
    
    for k = 1:nt
        q_dot = X(k,1:Nq).';
        v_dot = X(k,Nq+1:Nq+2).';
        q = X(k,Nq+3 : 2*Nq+2).';
        v = X(k,2*Nq+3:2*Nq+4).';

        Robot.setState(q, q_dot);
        Robot.Calculate_H_matrix();
        [B, D, K, H] = Robot.getDynamicMatrices();


        %Adding the end effector
        EEF = Robot.End_effector;

        [r_EEF,r_EEF_dot,~] = Robot.FW_Kinematics_2D(q, q_dot, 2, Robot.L(2));

        K_EEF = [EEF.Kx 0; 0 EEF.Ky];
        C_EEF = [EEF.Cx 0; 0 EEF.Cy];

        J = Robot.shape_functions();

        K_BASE = -J'*K_EEF;
        C_BASE = -J'*C_EEF;

        r_res    = r_EEF     - J*q;
        rdot_res = r_EEF_dot - J*q_dot;

        F_EEF  = -K_EEF*r_res - C_EEF*rdot_res;
        F_BASE =  J'*(K_EEF*r_res + C_EEF*rdot_res);

        K_ = [K+J'*K_EEF*J K_BASE;K_BASE' K_EEF];
        C_ = [D+J'*C_EEF*J C_BASE;C_BASE' C_EEF];
        M_ = blkdiag(B, EEF.M_pe * eye(2));
        H_ = [H;0;0] + [F_BASE;F_EEF];

        % Controller and current augmented matrices
        z     = [q; v];
        z_dot = [q_dot; v_dot];
        tau = compute_control_torque( ...
            t(k), z, z_dot, Robot, ...
            interp_q, interp_omega, interp_alpha, ...
            M_, C_, K_, H_);
        tau_history(k,:) = tau.';
        power_history(k,:) = (tau.*q_dot(1:Robot.N)).';
    end

end

function plot_state_histories(t, X, Robot)
    
    Nq = Robot.NDoFs;
    coordinate_names = strings(Nq+2,1);
    coordinate_names(end-1) = 'Vx';
    coordinate_names(end)   = 'Vy';
    
    for i = 1:Robot.N
        coordinate_names(i) = sprintf('theta_%d',i);
    end
    
    for i = 1:Robot.N
        for j = 1:Robot.m_i(i)
            idx = Robot.N + sum(Robot.m_i(1:i-1)) + j;
            coordinate_names(idx) = sprintf('delta_%d%d',i,j);
        end
    end
            coordinate_names(end) = sprintf('Y_bar');
    
    for i = 1:Nq+2
        figure('Color','w');
        plot(t,X(:,Nq+2+i),'LineWidth',1.6);
        grid on;
        xlabel('Time [s]');
        ylabel(coordinate_names(i));
        title(sprintf('%s versus time',coordinate_names(i)));
    end
    
    for i = 1:Nq+2
        figure('Color','w');
        plot(t,X(:,i),'LineWidth',1.6);
        grid on;
        xlabel('Time [s]');
        ylabel(sprintf('d%s/dt',coordinate_names(i)));
        title(sprintf('Velocity of %s',coordinate_names(i)));
    end

end

function plot_cartesian_path(Robot, actual_position, theta_ref)

    desired_position = rigid_reference_tip_positions(Robot,theta_ref);
    
    figure('Color','w','Name','Cartesian Path Tracking');
    plot(1000*desired_position(:,1),1000*desired_position(:,2), ...
        'r--','LineWidth',2,'DisplayName','Desired Path');
    hold on;
    plot(1000*actual_position(:,1),1000*actual_position(:,2), ...
        'b-', 'DisplayName','Actual Path');
    grid on;
    axis equal;
    xlabel('X [mm]');
    ylabel('Y [mm]');
    title('End-Effector Path Tracking');
    legend('Location','best');

end

function plot_tip_quantities(t, tip_position, tip_velocity, ...
                             EEF_position, EEF_velocity, ...
                             relative_position, relative_velocity)

    % Absolute position
    figure('Color','w','Name','End-Effector Position');

    tiledlayout(2,1);

    nexttile;
    plot(t,1000*tip_position(:,1),'LineWidth',1.5);
    hold on;
    plot(t,1000*EEF_position(:,1),'--','LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('X [mm]');
    legend('Robot tip r_{3x}','EEF v_x','Location','best');
    title('X Position');

    nexttile;
    plot(t,1000*tip_position(:,2),'LineWidth',1.5);
    hold on;
    plot(t,1000*EEF_position(:,2),'--','LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('Y [mm]');
    legend('Robot tip r_{3y}','EEF v_y','Location','best');
    title('Y Position');


    % Absolute velocity
    figure('Color','w','Name','End-Effector Velocity');

    tiledlayout(2,1);

    nexttile;
    plot(t,1000*tip_velocity(:,1),'LineWidth',1.5);
    hold on;
    plot(t,1000*EEF_velocity(:,1),'--','LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('V_x [mm/s]');
    legend('Robot tip','EEF','Location','best');
    title('X Velocity');

    nexttile;
    plot(t,1000*tip_velocity(:,2),'LineWidth',1.5);
    hold on;
    plot(t,1000*EEF_velocity(:,2),'--','LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('V_y [mm/s]');
    legend('Robot tip','EEF','Location','best');
    title('Y Velocity');


    % Relative position
    figure('Color','w','Name','Relative End-Effector Position');

    tiledlayout(2,1);

    nexttile;
    plot(t,1000*relative_position(:,1),'LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('v_x-r_{3x} [mm]');
    title('Relative X Position');

    nexttile;
    plot(t,1000*relative_position(:,2),'LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('v_y-r_{3y} [mm]');
    title('Relative Y Position');


    % Relative velocity
    figure('Color','w','Name','Relative End-Effector Velocity');

    tiledlayout(2,1);

    nexttile;
    plot(t,1000*relative_velocity(:,1),'LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('\dot{v}_x-\dot{r}_{3x} [mm/s]');
    title('Relative X Velocity');

    nexttile;
    plot(t,1000*relative_velocity(:,2),'LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('\dot{v}_y-\dot{r}_{3y} [mm/s]');
    title('Relative Y Velocity');

end

function plot_torque_power(t, tau, power)
    
    figure('Color','w');
    plot(t,tau,'LineWidth',1.6);
    grid on;
    xlabel('Time [s]');
    ylabel('Torque [Nm]');
    legend('Joint 1','Joint 2','Location','best');
    title('Joint Torques');
    
    figure('Color','w');
    plot(t,power,'LineWidth',1.6);
    grid on;
    xlabel('Time [s]');
    ylabel('Power [W]');
    legend('Joint 1','Joint 2','Location','best');
    title('Joint Power');

end

function desired_position = rigid_reference_tip_positions(Robot,theta_ref)
    
    n = size(theta_ref,2);
    desired_position = zeros(n,2);
    q = zeros(Robot.NDoFs,1);
    q_dot = zeros(Robot.NDoFs,1);
    
    for k = 1:n
        q(1:Robot.N) = theta_ref(:,k);
        [desired_position(k,:),~] = Robot.FW_Kinematics_2D( ...
            q, q_dot, Robot.N, Robot.L(Robot.N));
    end

end

function live = init_live_plots(Robot,theta_ref,t_span)
    
    desired_position = rigid_reference_tip_positions(Robot,theta_ref);
    live = struct();
    
    live.figPath = figure('Name','Live Cartesian Path','Color','w');
    live.axPath = axes('Parent',live.figPath);
    hold(live.axPath,'on');
    grid(live.axPath,'on');
    axis(live.axPath,'equal');
    plot(live.axPath,1000*desired_position(:,1),1000*desired_position(:,2), ...
        'r--','LineWidth',2);
    live.pathActual = animatedline(live.axPath,'Color','Blue');
    live.pathTip = plot(live.axPath,nan,nan,'ko','MarkerFaceColor','k');
    xlabel(live.axPath,'X [mm]');
    ylabel(live.axPath,'Y [mm]');
    legend(live.axPath,{'Desired path','Actual path','Current tip'}, ...
        'Location','best');
    live.statusText = text(live.axPath,0.02,0.98,'Preparing solver ...', ...
        'Units','normalized','VerticalAlignment','top', ...
        'BackgroundColor','w','Margin',5);
    live.tStart = t_span(1);
    live.tEnd = t_span(end);

end

function status = live_plot_callback(t,x,flag,Robot,live)
    
    status = 0;
    persistent wallTimer lastPrintedPercent outputPointCount
    
    if strcmp(flag,'init')
        wallTimer = tic;
        lastPrintedPercent = 0;
        outputPointCount = 0;
        set(live.statusText,'String','Starting integration: 0.0%');
        fprintf('Simulation started. Progress will be reported every 5%%.\n');
        drawnow;
        return;
    elseif strcmp(flag,'done')
        elapsed = toc(wallTimer);
        set(live.statusText,'String',sprintf( ...
            'Integration complete | wall time %.1f s',elapsed));
        fprintf('Simulation integration complete in %.1f s.\n',elapsed);
        drawnow;
        return;
    elseif ~isempty(flag)
        return;
    end
    
    if ~isgraphics(live.figPath)
        status = 1;
        fprintf('Simulation stopped because the live-path window was closed.\n');
        return;
    end
    
    if any(~isfinite(x(:)))
        warning('LightWeightRoboticArmSims:NonFiniteState', ...
            'The solver produced NaN or Inf. Integration was stopped.');
        status = 1;
        return;
    end
    
    Nq = Robot.NDoFs;
    outputPointCount = outputPointCount + numel(t);
    
    for k = 1:numel(t)
        q_dot = x(1:Nq,k);
        q = x(Nq+3:2*Nq+2,k);
    
        [position,~] = Robot.FW_Kinematics_2D( ...
            q, q_dot, Robot.N, Robot.L(Robot.N));
    
        addpoints(live.pathActual,1000*position(1),1000*position(2));
        set(live.pathTip,'XData',1000*position(1),'YData',1000*position(2));
    end
    
    progress = 100*(t(end)-live.tStart)/(live.tEnd-live.tStart);
    progress = max(0,min(100,progress));
    elapsed = toc(wallTimer);
    if progress >= 0.5
        remaining = elapsed*(100-progress)/progress;
        etaText = sprintf('%.1f s',remaining);
    else
        remaining = NaN;
        etaText = 'estimating';
    end
    
    qCurrent = x(Nq+3:2*Nq+2,end);
    qDotCurrent = x(1:Nq,end);
    if Nq > Robot.N
        maxFlexible = max(abs(qCurrent(Robot.N+1:end)));
    else
        maxFlexible = 0;
    end
    
    set(live.statusText,'String',sprintf([ ...
        't = %.4f / %.4f s | %.1f%%\n' ...
        'wall %.1f s | ETA %s | output points %d\n' ...
        'max |flexible q| %.3e | max |q dot| %.3e'], ...
        t(end),live.tEnd,progress,elapsed,etaText,outputPointCount, ...
        maxFlexible,max(abs(qDotCurrent))));
    
    reportedPercent = 5*floor(progress/5);
    if reportedPercent >= lastPrintedPercent + 5
        fprintf('Progress %3.0f%%: simulated t = %.4f s, wall = %.1f s, ETA = %.1f s.\n', ...
            reportedPercent,t(end),elapsed,remaining);
        lastPrintedPercent = reportedPercent;
    end
    
    drawnow limitrate;

end

function animate_robot(Robot,t,X,theta_ref)
    
    desired_position = rigid_reference_tip_positions(Robot,theta_ref);
    reach = sum(Robot.L);
    
    figure('Color','w','Name','Robot Animation with Path Tracking');
    plot(desired_position(:,1),desired_position(:,2),'r--','LineWidth',1.5);
    hold on;
    grid on;
    axis equal;
    axis(1.05*[-reach reach -reach reach]);
    
    Nq = Robot.NDoFs;
    frame_step = max(1,round(numel(t)/250));
    h_robot = gobjects(1);
    h_joints = gobjects(1);
    
    for k = 1:frame_step:numel(t)
        if isgraphics(h_robot)
            delete(h_robot);
        end
    
        if isgraphics(h_joints)
            delete(h_joints);
        end
    
        q_dot = X(k,1:Nq).';
        q = X(k,Nq+3:2*Nq+2).';
    
        nodes = zeros(Robot.N+1,2);
        for i = 1:Robot.N
            [nodes(i+1,:),~] = Robot.FW_Kinematics_2D(q,q_dot,i,Robot.L(i));
        end
    
        h_robot = plot(nodes(:,1),nodes(:,2),'b-','LineWidth',2);
        h_joints = plot(nodes(:,1),nodes(:,2),'ko','MarkerFaceColor','k');
        title(sprintf('Time: %.3f s',t(k)));
        drawnow limitrate;
    end

end


function plot_relative_stft(t, relative_position, freq_mean, freq_sd)

% ============================================================
% Relative-position frequency analysis
%
% INPUTS:
%   t                 - simulation time vector
%   relative_position - [x_rel, y_rel]
%   freq_mean         - mean expected modal frequencies [Hz]
%   freq_sd           - standard deviation of modal frequencies [Hz]
%
% Produces:
%   1. STFT of X relative displacement [dB]
%   2. STFT of Y relative displacement [dB]
%   3. Combined X-Y STFT [dB]
%   4. Full-record amplitude spectrum
%   5. Full-record spectrum in dB
%
% Expected modal frequencies are added to the frequency axes
% as:
%
%        mean +/- standard deviation
%
% ============================================================


% ------------------------------------------------------------
% Optional frequency inputs
% ------------------------------------------------------------

if nargin < 3
    freq_mean = [];
end

if nargin < 4
    freq_sd = zeros(size(freq_mean));
end


freq_mean = freq_mean(:).';
freq_sd   = freq_sd(:).';


% Mean and standard deviation vectors must match
if length(freq_mean) ~= length(freq_sd)

    error(['freq_mean and freq_sd must contain the same ' ...
           'number of elements.']);

end


% ------------------------------------------------------------
% Remove invalid modal-frequency values
% ------------------------------------------------------------

valid = ...
    isfinite(freq_mean) & ...
    isfinite(freq_sd) & ...
    freq_mean >= 0 & ...
    freq_sd >= 0;

freq_mean = freq_mean(valid);
freq_sd   = freq_sd(valid);


% ------------------------------------------------------------
% Sort modes according to mean frequency
% ------------------------------------------------------------

[freq_mean,idx] = sort(freq_mean);

freq_sd = freq_sd(idx);



% ============================================================
% INPUT FORMATTING
% ============================================================

t = t(:);

ex_raw = relative_position(:,1);
ey_raw = relative_position(:,2);



% ============================================================
% UNIFORM RESAMPLING
% ============================================================

dt_solver = median(diff(t));

% Do not artificially upsample solver data.
%
% If solver sampling is faster than 1 ms:
%     use 1 ms.
%
% If solver sampling is slower:
%     keep the solver sampling period.

dt = max(dt_solver,1e-3);

Fs = 1/dt;


t_uniform = (t(1):dt:t(end)).';


ex = interp1( ...
    t, ...
    ex_raw, ...
    t_uniform, ...
    'pchip');

ey = interp1( ...
    t, ...
    ey_raw, ...
    t_uniform, ...
    'pchip');



% ============================================================
% REMOVE MEAN / DC COMPONENT
% ============================================================

ex = ex - mean(ex);

ey = ey - mean(ey);



% ============================================================
% STFT SETTINGS
% ============================================================

totalTime = t_uniform(end) - t_uniform(1);


% ------------------------------------------------------------
% Longer window improves modal-frequency resolution
% ------------------------------------------------------------

windowDuration = ...
    min(2.0,totalTime/2);


windowLength = ...
    round(windowDuration*Fs);


% Minimum allowed window
windowLength = ...
    max(32,windowLength);


% Window cannot exceed signal length
windowLength = ...
    min(windowLength,length(ex));


window = ...
    hann(windowLength,'periodic');


% 75 percent overlap
overlap = ...
    round(0.75*windowLength);


overlap = ...
    min(overlap,windowLength-1);


% ------------------------------------------------------------
% FFT length
%
% Zero padding improves graphical frequency sampling,
% but does NOT improve true frequency resolution.
% ------------------------------------------------------------

nfft = ...
    2^nextpow2(4*windowLength);


% Approximate actual frequency resolution
frequencyResolution = ...
    Fs/windowLength;



% ============================================================
% DISPLAY ANALYSIS INFORMATION
% ============================================================

fprintf('\n');

fprintf('=====================================================\n');
fprintf(' Relative-position spectral analysis\n');
fprintf('=====================================================\n');

fprintf('Simulation duration      : %.4f s\n', ...
    totalTime);

fprintf('Uniform sample time      : %.6f s\n', ...
    dt);

fprintf('Sampling frequency       : %.2f Hz\n', ...
    Fs);

fprintf('STFT window duration     : %.4f s\n', ...
    windowLength/Fs);

fprintf('Approx. frequency res.   : %.4f Hz\n', ...
    frequencyResolution);

fprintf('FFT length               : %d\n', ...
    nfft);



% ------------------------------------------------------------
% Print expected modal frequencies
% ------------------------------------------------------------

if ~isempty(freq_mean)

    fprintf('\n');
    fprintf('Expected modal frequencies:\n');

    for k = 1:length(freq_mean)

        fprintf( ...
            '  f_%d = %.4f +/- %.4f Hz\n', ...
            k, ...
            freq_mean(k), ...
            freq_sd(k));

    end

end


fprintf('=====================================================\n');
fprintf('\n');



% ============================================================
% X RELATIVE DISPLACEMENT STFT
% ============================================================

[Sx,F,Tx] = spectrogram( ...
    ex, ...
    window, ...
    overlap, ...
    nfft, ...
    Fs);


Ax = abs(Sx)/sum(window);



% ------------------------------------------------------------
% Convert one-sided FFT magnitude to amplitude
% ------------------------------------------------------------

if mod(nfft,2) == 0

    Ax(2:end-1,:) = ...
        2*Ax(2:end-1,:);

else

    Ax(2:end,:) = ...
        2*Ax(2:end,:);

end



% ============================================================
% Y RELATIVE DISPLACEMENT STFT
% ============================================================

[Sy,~,Ty] = spectrogram( ...
    ey, ...
    window, ...
    overlap, ...
    nfft, ...
    Fs);


Ay = abs(Sy)/sum(window);



% ------------------------------------------------------------
% Convert one-sided FFT magnitude to amplitude
% ------------------------------------------------------------

if mod(nfft,2) == 0

    Ay(2:end-1,:) = ...
        2*Ay(2:end-1,:);

else

    Ay(2:end,:) = ...
        2*Ay(2:end,:);

end



% ============================================================
% COMBINED X-Y AMPLITUDE
% ============================================================

Arel = ...
    sqrt(Ax.^2 + Ay.^2);



% ============================================================
% CONVERT STFT AMPLITUDES TO dB
% ============================================================

Ax_dB = ...
    20*log10( ...
    Ax/(max(Ax(:)) + eps) + eps);


Ay_dB = ...
    20*log10( ...
    Ay/(max(Ay(:)) + eps) + eps);


Arel_dB = ...
    20*log10( ...
    Arel/(max(Arel(:)) + eps) + eps);



% ============================================================
% X STFT
% ============================================================

figure( ...
    'Color','w', ...
    'Name','STFT Relative X Position');


meshc( ...
    Tx, ...
    F, ...
    Ax_dB);


xlabel('Time [s]');

ylabel('Frequency [Hz]');

zlabel('Relative Amplitude [dB]');


title('STFT of Relative X Displacement');


grid on;

axis tight;


zlim([-80 0]);


view(45,35);


colorbar;



% ------------------------------------------------------------
% Add expected modal frequencies to Y-axis ticks
% ------------------------------------------------------------

if ~isempty(freq_mean)

    valid = ...
        freq_mean <= max(F);


    fm = freq_mean(valid);

    fs = freq_sd(valid);


    currentTicks = yticks;


    newTicks = sort(unique( ...
        [currentTicks fm]));


    yticks(newTicks);


    labels = strings(size(newTicks));


    % Standard MATLAB ticks
    for k = 1:length(newTicks)

        labels(k) = ...
            sprintf('%.2f',newTicks(k));

    end


    % Replace expected modal-frequency ticks
    for k = 1:length(fm)

        [~,ind] = ...
            min(abs(newTicks - fm(k)));


        labels(ind) = ...
            sprintf( ...
            '%.2f +/- %.2f', ...
            fm(k), ...
            fs(k));

    end


    yticklabels(labels);

end



% ============================================================
% Y STFT
% ============================================================

figure( ...
    'Color','w', ...
    'Name','STFT Relative Y Position');


meshc( ...
    Ty, ...
    F, ...
    Ay_dB);


xlabel('Time [s]');

ylabel('Frequency [Hz]');

zlabel('Relative Amplitude [dB]');


title('STFT of Relative Y Displacement');


grid on;

axis tight;


zlim([-80 0]);


view(45,35);


colorbar;



% ------------------------------------------------------------
% Add expected modal frequencies to Y-axis ticks
% ------------------------------------------------------------

if ~isempty(freq_mean)

    valid = ...
        freq_mean <= max(F);


    fm = freq_mean(valid);

    fs = freq_sd(valid);


    currentTicks = yticks;


    newTicks = sort(unique( ...
        [currentTicks fm]));


    yticks(newTicks);


    labels = strings(size(newTicks));


    for k = 1:length(newTicks)

        labels(k) = ...
            sprintf('%.2f',newTicks(k));

    end


    for k = 1:length(fm)

        [~,ind] = ...
            min(abs(newTicks - fm(k)));


        labels(ind) = ...
            sprintf( ...
            '%.2f +/- %.2f', ...
            fm(k), ...
            fs(k));

    end


    yticklabels(labels);

end



% ============================================================
% COMBINED X-Y STFT
% ============================================================

figure( ...
    'Color','w', ...
    'Name','STFT Combined Relative Position');


meshc( ...
    Tx, ...
    F, ...
    Arel_dB);


xlabel('Time [s]');

ylabel('Frequency [Hz]');

zlabel('Relative Amplitude [dB]');


title('Combined X-Y Relative Position STFT');


grid on;

axis tight;


zlim([-80 0]);


view(45,35);


colorbar;



% ------------------------------------------------------------
% Add expected modal frequencies to Y-axis ticks
% ------------------------------------------------------------

if ~isempty(freq_mean)

    valid = ...
        freq_mean <= max(F);


    fm = freq_mean(valid);

    fs = freq_sd(valid);


    currentTicks = yticks;


    newTicks = sort(unique( ...
        [currentTicks fm]));


    yticks(newTicks);


    labels = strings(size(newTicks));


    for k = 1:length(newTicks)

        labels(k) = ...
            sprintf('%.2f',newTicks(k));

    end


    for k = 1:length(fm)

        [~,ind] = ...
            min(abs(newTicks - fm(k)));


        labels(ind) = ...
            sprintf( ...
            '%.2f +/- %.2f', ...
            fm(k), ...
            fs(k));

    end


    yticklabels(labels);

end



% ============================================================
% FULL-RECORD FFT
% ============================================================

N = length(t_uniform);


Nfft_full = ...
    2^nextpow2(N);


X = fft( ...
    ex, ...
    Nfft_full);


Y = fft( ...
    ey, ...
    Nfft_full);


f_full = ...
    Fs*(0:(Nfft_full/2))/Nfft_full;



% ============================================================
% X AMPLITUDE SPECTRUM
% ============================================================

Px = abs(X/N);


Px = ...
    Px(1:Nfft_full/2+1);


if length(Px) > 2

    Px(2:end-1) = ...
        2*Px(2:end-1);

end



% ============================================================
% Y AMPLITUDE SPECTRUM
% ============================================================

Py = abs(Y/N);


Py = ...
    Py(1:Nfft_full/2+1);


if length(Py) > 2

    Py(2:end-1) = ...
        2*Py(2:end-1);

end



% ============================================================
% COMBINED X-Y SPECTRUM
% ============================================================

Pxy = ...
    sqrt(Px.^2 + Py.^2);



% ============================================================
% FULL-RECORD AMPLITUDE SPECTRUM
% ============================================================

figure( ...
    'Color','w', ...
    'Name','Relative Position Full Spectrum');


plot( ...
    f_full, ...
    1000*Px, ...
    'LineWidth',1.2);


hold on;


plot( ...
    f_full, ...
    1000*Py, ...
    'LineWidth',1.2);


plot( ...
    f_full, ...
    1000*Pxy, ...
    'LineWidth',1.5);


hold off;


xlabel('Frequency [Hz]');

ylabel('Amplitude [mm]');


title('Full-Record Relative-Position Spectrum');


legend( ...
    'X relative', ...
    'Y relative', ...
    'Combined X-Y', ...
    'Location','best');


grid on;

axis tight;



% ------------------------------------------------------------
% Add expected modal frequencies to X-axis ticks
% ------------------------------------------------------------

if ~isempty(freq_mean)

    valid = ...
        freq_mean <= max(f_full);


    fm = freq_mean(valid);

    fs = freq_sd(valid);


    currentTicks = xticks;


    newTicks = sort(unique( ...
        [currentTicks fm]));


    xticks(newTicks);


    labels = strings(size(newTicks));


    for k = 1:length(newTicks)

        labels(k) = ...
            sprintf('%.2f',newTicks(k));

    end


    for k = 1:length(fm)

        [~,ind] = ...
            min(abs(newTicks - fm(k)));


        labels(ind) = ...
            sprintf( ...
            '%.2f +/- %.2f', ...
            fm(k), ...
            fs(k));

    end


    xticklabels(labels);

end



% ============================================================
% FULL-RECORD SPECTRUM IN dB
% ============================================================

Pxy_dB = ...
    20*log10( ...
    Pxy/(max(Pxy) + eps) + eps);



figure( ...
    'Color','w', ...
    'Name','Relative Position Spectrum dB');


plot( ...
    f_full, ...
    Pxy_dB, ...
    'LineWidth',1.5);


xlabel('Frequency [Hz]');

ylabel('Relative Amplitude [dB]');


title('Combined Relative-Position Spectrum');


grid on;

axis tight;


ylim([-100 0]);



% ------------------------------------------------------------
% Add expected modal frequencies to X-axis ticks
% ------------------------------------------------------------

if ~isempty(freq_mean)

    valid = ...
        freq_mean <= max(f_full);


    fm = freq_mean(valid);

    fs = freq_sd(valid);


    currentTicks = xticks;


    newTicks = sort(unique( ...
        [currentTicks fm]));


    xticks(newTicks);


    labels = strings(size(newTicks));


    for k = 1:length(newTicks)

        labels(k) = ...
            sprintf('%.2f',newTicks(k));

    end


    for k = 1:length(fm)

        [~,ind] = ...
            min(abs(newTicks - fm(k)));


        labels(ind) = ...
            sprintf( ...
            '%.2f +/- %.2f', ...
            fm(k), ...
            fs(k));

    end


    xticklabels(labels);

end


end

function plot_absolute_stft(t, absolute_position)

% ------------------------------------------------------------
% Uniform resampling
% ------------------------------------------------------------

t = t(:);

dt_solver = median(diff(t));
dt = min(dt_solver,1e-3);

Fs = 1/dt;

t_uniform = (t(1):dt:t(end)).';

x_abs = interp1(t,absolute_position(:,1), ...
                t_uniform,'pchip');

y_abs = interp1(t,absolute_position(:,2), ...
                t_uniform,'pchip');


% ------------------------------------------------------------
% Remove mean / DC component
% ------------------------------------------------------------

x_abs = x_abs - mean(x_abs);
y_abs = y_abs - mean(y_abs);


% ------------------------------------------------------------
% STFT settings
% ------------------------------------------------------------

totalTime = t_uniform(end) - t_uniform(1);

windowDuration = min(0.25,totalTime/4);

windowLength = max(32,round(windowDuration*Fs));

window = hann(windowLength,'periodic');

overlap = round(0.75*windowLength);

nfft = 2^nextpow2(4*windowLength);


% ------------------------------------------------------------
% X absolute displacement
% ------------------------------------------------------------

[Sx,F,Tx] = spectrogram( ...
    x_abs,window,overlap,nfft,Fs);

Ax = abs(Sx)/sum(window);

if mod(nfft,2) == 0
    Ax(2:end-1,:) = 2*Ax(2:end-1,:);
else
    Ax(2:end,:) = 2*Ax(2:end,:);
end


% ------------------------------------------------------------
% Y absolute displacement
% ------------------------------------------------------------

[Sy,F,Ty] = spectrogram( ...
    y_abs,window,overlap,nfft,Fs);

Ay = abs(Sy)/sum(window);

if mod(nfft,2) == 0
    Ay(2:end-1,:) = 2*Ay(2:end-1,:);
else
    Ay(2:end,:) = 2*Ay(2:end,:);
end


% ------------------------------------------------------------
% Plot X absolute-position STFT
% ------------------------------------------------------------

figure('Color','w','Name','STFT Absolute X Position');

meshc(Tx,F,1000*Ax);

xlabel('Time [s]');
ylabel('Frequency [Hz]');
zlabel('Amplitude [mm]');

title('STFT of Absolute EEF X Position');

grid on;
axis tight;

view(45,35);

colorbar;


% ------------------------------------------------------------
% Plot Y absolute-position STFT
% ------------------------------------------------------------

figure('Color','w','Name','STFT Absolute Y Position');

meshc(Ty,F,1000*Ay);

xlabel('Time [s]');
ylabel('Frequency [Hz]');
zlabel('Amplitude [mm]');

title('STFT of Absolute EEF Y Position');

grid on;
axis tight;

view(45,35);

colorbar;

end