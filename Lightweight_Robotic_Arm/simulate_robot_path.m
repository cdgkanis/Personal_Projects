function [theta, omega, alpha] = simulate_robot_path(Robot, tau, F, mov_params)
    % simulate_robot_path Calculates inverse kinematics, joint velocities, and accelerations.
    
    % Scaling factor for acceleration (since time is squared in the denominator)
    F_dot = F/(tau(end) - tau(1));

    % ==========================================
    % 1. Arm Parameters
    % ==========================================
    L1 = Robot.L(1);
    L2 = Robot.L(2);
    R_robot_min = 0.2; % Minimum reach constraint

    % ==========================================
    % 2. Generate Target Trajectory
    % ==========================================
    [X_val, Y_val] = Path(tau, mov_params);
    target_x = X_val;
    target_y = Y_val;
    
    % ==========================================
    % 3. Validate Workspace
    % ==========================================
    [valid_mask, status] = validate_workspace(target_x, target_y, L1, L2, R_robot_min);
    
    % Separate the path into reachable and unreachable points
    valid_x = target_x(valid_mask);
    valid_y = target_y(valid_mask);
    
    % Preallocate outputs in case the entire path is invalid
    theta = [];
    omega = [];
    alpha = [];
    
    % ==========================================
    % 4. Inverse Kinematics (Positions)
    % ==========================================
    if any(valid_mask)
        [theta1_rad, theta2_rad] = calculate_inverse_kinematics(valid_x, valid_y, L1, L2);
        
        % Assign theta output: 2xN array (Row 1: Base, Row 2: Elbow)
        theta = [theta1_rad; theta2_rad];
    end
    
    % ==========================================
    % 5. Exact Joint Velocities & Accelerations
    % ==========================================
    if any(valid_mask)
        valid_tau = tau(valid_mask);
        
        [X_dot_val, Y_dot_val] = Path_dot(valid_tau, mov_params);
        [X_ddot_val, Y_ddot_val] = Path_ddot(valid_tau, mov_params);
        
        % Scale to real-world Cartesian velocities and accelerations
        valid_x_dot = F * X_dot_val;
        valid_y_dot = F * Y_dot_val;
        
        valid_x_ddot = F_dot * X_ddot_val;
        valid_y_ddot = F_dot * Y_ddot_val;
        
        % Preallocate velocity and acceleration arrays
        omega1 = zeros(1, length(valid_x));
        omega2 = zeros(1, length(valid_x));
        alpha1 = zeros(1, length(valid_x));
        alpha2 = zeros(1, length(valid_x));
        
        for i = 1:length(valid_x)
            th1 = theta1_rad(i);
            th2 = theta2_rad(i);
            
            % --- A. VELOCITY ---
            J = [-L1*sin(th1) - L2*sin(th1+th2),  -L2*sin(th1+th2);
                  L1*cos(th1) + L2*cos(th1+th2),   L2*cos(th1+th2)];
              
            V_cartesian = [valid_x_dot(i); valid_y_dot(i)];
            V_joint = J \ V_cartesian;
            
            om1 = V_joint(1);
            om2 = V_joint(2);
            omega1(i) = om1;
            omega2(i) = om2;
            
            % --- B. ACCELERATION ---
            J11_dot = -L1*cos(th1)*om1 - L2*cos(th1+th2)*(om1+om2);
            J12_dot = -L2*cos(th1+th2)*(om1+om2);
            J21_dot = -L1*sin(th1)*om1 - L2*sin(th1+th2)*(om1+om2);
            J22_dot = -L2*sin(th1+th2)*(om1+om2);
            
            J_dot = [J11_dot, J12_dot;
                     J21_dot, J22_dot];
                 
            A_cartesian = [valid_x_ddot(i); valid_y_ddot(i)];
            A_joint = J \ (A_cartesian - J_dot * V_joint);
            
            alpha1(i) = A_joint(1);
            alpha2(i) = A_joint(2);
        end
        
        % Assign outputs: 2xN arrays
        omega = [omega1; omega2];
        alpha = [alpha1; alpha2];
    end
    
    % ==========================================
    % 6. Call Extracted Plotting/Printing Function
    % ==========================================
    plot_robot_results(tau, valid_mask, status, omega, alpha, target_x, target_y, Robot, F, F_dot);
end

function plot_robot_results(tau, valid_mask, status, omega, alpha, target_x, target_y, Robot, F, F_dot)
    % PLOT_ROBOT_RESULTS Handles all text display and graphical plotting 
    % for the 2-DOF robotic arm path simulation.
    
    % ==========================================
    % 1. Print Workspace Validation Status
    % ==========================================
    disp(status);
    
    % Extract parameters and slice arrays
    L1 = Robot.L(1);
    L2 = Robot.L(2);
    R_robot_min = 0.2; % Minimum reach constraint
    
    valid_tau = tau(valid_mask);
    valid_x = target_x(valid_mask);
    valid_y = target_y(valid_mask);
    invalid_x = target_x(~valid_mask);
    invalid_y = target_y(~valid_mask);

    % ==========================================
    % 2. Plot Kinematics (Velocities & Accelerations)
    % ==========================================
    if any(valid_mask) && ~isempty(omega) && ~isempty(alpha)
        omega1 = omega(1, :);
        omega2 = omega(2, :);
        alpha1 = alpha(1, :);
        alpha2 = alpha(2, :);
        
        % --- PLOT VELOCITIES ---
        figure('Name', 'Joint Velocities', 'Color', 'w');
        plot(valid_tau, omega1, 'b', 'LineWidth', 1.5, 'DisplayName', '\omega_1 (Base)');
        hold on; grid on;
        plot(valid_tau, omega2, 'r', 'LineWidth', 1.5, 'DisplayName', '\omega_2 (Elbow)');
        title(sprintf('Exact Joint Velocities (F^*=%.2f, F\\_dot^*=%.2f)', F, F_dot));
        xlabel('Normalized Time (\tau)'); ylabel('Angular Velocity (rad/s)');
        legend('Location', 'best');
        
        % --- PLOT ACCELERATIONS ---
        figure('Name', 'Joint Accelerations', 'Color', 'w');
        plot(valid_tau, alpha1, 'b', 'LineWidth', 1.5, 'DisplayName', '\alpha_1 (Base)');
        hold on; grid on;
        plot(valid_tau, alpha2, 'r', 'LineWidth', 1.5, 'DisplayName', '\alpha_2 (Elbow)');
        title(sprintf('Exact Joint Accelerations (F^*=%.2f, F\\_dot^*=%.2f)', F, F_dot));
        xlabel('Normalized Time (\tau)'); ylabel('Angular Accel (rad/s^2)');
        legend('Location', 'best');
    end
    
    % ==========================================
    % 3. Plot Workspace and Path Boundaries
    % ==========================================
    figure('Name', 'Workspace Validation', 'Color', 'w');
    hold on; grid on; axis equal;
    title(sprintf('Workspace Validation (F^*=%.2f, F\\_dot^*=%.2f)', F, F_dot));
    xlabel('X Coordinate (m)'); ylabel('Y Coordinate (m)');
    
    circle_theta = linspace(0, 2*pi, 100);
    
    % Draw Outer Boundary
    plot((L1+L2)*cos(circle_theta), (L1+L2)*sin(circle_theta), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Max Reach (R_{max})');
    
    % Draw Inner Boundary
    inner_radius = max([abs(L1-L2), R_robot_min]);
    if inner_radius > 0
        plot(inner_radius * cos(circle_theta), inner_radius * sin(circle_theta), 'k:', 'LineWidth', 1.5, 'DisplayName', 'Min Reach (R_{min})');
    end
    
    % Scatter points
    if ~isempty(valid_x)
        scatter(valid_x, valid_y, 20, 'g', 'filled', 'DisplayName', 'Reachable Points');
    end
    if ~isempty(invalid_x)
        scatter(invalid_x, invalid_y, 20, 'r', 'filled', 'DisplayName', 'Unreachable Points');
    end
    
    scatter(0, 0, 100, 'k', 'filled', 'DisplayName', 'Base Origin');
    legend('Location', 'best');
end