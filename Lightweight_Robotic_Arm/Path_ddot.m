function [x_ddot, y_ddot] = Path_ddot(t, mov_params)
    % Straight edges => zero acceleration except impulsive corner transitions,
    % which are not represented in continuous-time output.
    x_ddot = zeros(size(t));
    y_ddot = zeros(size(t));
end