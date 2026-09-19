function [x_dot, y_dot] = Path_dot(t, mov_params)
%PATH_DOT Edge velocity of the unsmoothed path returned by Path.
%
% Velocity changes discontinuously at the three geometric corners. At an
% exact corner sample, this function returns the direction of the segment
% selected by the same inequalities used in Path.

    [A, B, C, dAB, dBC, dCA, P, s1, s2] = geometry(mov_params);
    [s, sdot] = pathTiming(t, mov_params.T, P, isBidirectional(mov_params));

    x_dot = zeros(size(t));
    y_dot = zeros(size(t));

    uAB = (B - A)/dAB;
    uBC = (C - B)/dBC;
    uCA = (A - C)/dCA;

    for i = 1:numel(t)
        si = s(i);

        if si < s1
            direction = uAB;
        elseif si < s2
            direction = uBC;
        else
            direction = uCA;
        end

        velocity = direction*sdot(i);
        x_dot(i) = velocity(1);
        y_dot(i) = velocity(2);
    end
end

function [s, sdot] = pathTiming(t, T, P, bidirectional)
    if ~(isscalar(T) && isfinite(T) && T > 0)
        error('Path_dot:InvalidPeriod', ...
              'mov_params.T must be a positive scalar.');
    end

    tau = mod(t, T);
    if bidirectional
        phase = 2*pi*tau/T;
        frac = 0.5*(1 - cos(phase));
        fracDot = (pi/T)*sin(phase);
    else
        frac = tau/T;
        fracDot = ones(size(t))/T;
    end

    s = P*frac;
    sdot = P*fracDot;
end

function value = isBidirectional(mov_params)
    value = isfield(mov_params, 'bidirectional') && ...
            logical(mov_params.bidirectional);
end

function [A, B, C, dAB, dBC, dCA, P, s1, s2] = geometry(mov_params)
    a = mov_params.a;
    if ~(isscalar(a) && isfinite(a) && a > 0)
        error('Path_dot:InvalidScale', ...
              'mov_params.a must be a positive scalar.');
    end

    A = [0, 0];
    B = [a, 0];
    C = [0, a*sqrt(3)];
    dAB = norm(B - A);
    dBC = norm(C - B);
    dCA = norm(A - C);
    P = dAB + dBC + dCA;
    s1 = dAB;
    s2 = dAB + dBC;
end
