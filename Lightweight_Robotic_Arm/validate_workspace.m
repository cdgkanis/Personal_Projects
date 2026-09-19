function [is_reachable, status_msg] = validate_workspace(x, y, L1, L2, R_robot_min)
    % Calculates the radial distance of each point from the base origin
    r = sqrt(x.^2 + y.^2);
    
    % Define the physical boundaries of the arm
    R_max = L1 + L2;
    R_min = max([abs(L1 - L2), R_robot_min]);
    
    % Check if each point falls within the allowable limits
    is_reachable = (r >= R_min) & (r <= R_max);
    
    % Generate a summary status message
    if all(is_reachable)
        status_msg = 'Success: The entire path is inside the workspace.';
    elseif any(is_reachable)
        status_msg = 'Warning: A portion of the path is outside the workspace!';
    else
        status_msg = 'Error: The entire path is outside the workspace!';
    end
end