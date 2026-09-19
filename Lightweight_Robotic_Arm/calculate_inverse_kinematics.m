function [theta1, theta2] = calculate_inverse_kinematics(x, y, L1, L2)
    % Calculates joint angles (theta1, theta2) for a 2-DOF planar robotic arm
    % Inputs:
    %   x, y         - Arrays of target Cartesian coordinates
    %   L1, L2       - Link lengths
    %   elbow_config - String: 'down' or 'up'
    % Outputs:
    %   theta1, theta2 - Arrays of joint angles in radians

    % Calculate the squared distance from base to target
    r_sq = x.^2 + y.^2;
    
    % Law of Cosines to find the angle of the second joint
    D = (r_sq - L1^2 - L2^2) ./ (2 * L1 * L2);
    
    % Clip D to [-1, 1] to prevent complex numbers if a point is slightly out of reach
    D = max(min(D, 1.0), -1.0);
    
    % Calculate theta2 based on elbow configuration
    theta2 = atan2(sqrt(1 - D.^2), D);
    
    % Calculate theta1 using the result of theta2
    theta1 = atan2(y, x) - atan2(L2 .* sin(theta2), L1 + L2 .* cos(theta2));
    
    % Unwrap angles to prevent 2*pi jumps which cause massive torque spikes in PD controllers
    theta1 = unwrap(theta1);
    theta2 = unwrap(theta2);
end