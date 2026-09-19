function [x, y] = Path(t, mov_params)
%PATH Piecewise-linear triangular path without corner smoothing.
%
% Required mov_params fields:
%   a, cx, cy, T
%
% Optional mov_params field:
%   bidirectional : false (default) gives one traversal per period.
%                   true gives one forward and one reverse traversal per
%                   period, with zero speed at both reversal points.
%
% corner_radius and sigmoid_k are accepted in mov_params for interface
% compatibility but are not used by this unsmoothed path.

    [A, B, C, dAB, dBC, dCA, P, s1, s2] = geometry(mov_params);
    s = pathPosition(t, mov_params.T, P, isBidirectional(mov_params));

    x = zeros(size(t));
    y = zeros(size(t));

    uAB = (B - A)/dAB;
    uBC = (C - B)/dBC;
    uCA = (A - C)/dCA;

    for i = 1:numel(t)
        si = s(i);

        if si < s1
            pt = A + si*uAB;
        elseif si < s2
            pt = B + (si - s1)*uBC;
        else
            pt = C + (si - s2)*uCA;
        end

        x(i) = pt(1) + mov_params.cx;
        y(i) = pt(2) + mov_params.cy;
    end
end

function s = pathPosition(t, T, P, bidirectional)
    if ~(isscalar(T) && isfinite(T) && T > 0)
        error('Path:InvalidPeriod', 'mov_params.T must be a positive scalar.');
    end

    tau = mod(t, T);
    if bidirectional
        frac = 0.5*(1 - cos(2*pi*tau/T));
    else
        frac = tau/T;
    end
    s = P*frac;
end

function value = isBidirectional(mov_params)
    value = isfield(mov_params, 'bidirectional') && ...
            logical(mov_params.bidirectional);
end

function [A, B, C, dAB, dBC, dCA, P, s1, s2] = geometry(mov_params)
    a = mov_params.a;
    if ~(isscalar(a) && isfinite(a) && a > 0)
        error('Path:InvalidScale', 'mov_params.a must be a positive scalar.');
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
