function [Tau1, Tau2] = LQR(M, D, K, H, z, z_dot, ...
        q_d, q_dot_d, q_ddot_d, Q, R)
%T Reduced-order LQR trajectory controller for two rigid joints.
%
% LQR state:
%   x_e = [q_r - q_r,d; q_dot_r - q_dot_r,d]
%
% z = [q_robot; v_EEF] is the complete augmented coordinate vector.
% Only the first two rigid coordinates enter the feedback error; all other
% coordinates remain in the inverse-dynamics feedforward calculation.
%
% Q must be 4-by-4 and R must be 2-by-2.  The flexible coordinates are
% excluded from the feedback state, but their dynamic coupling is retained
% in the Schur-complement inverse-dynamics feedforward term.

    nr = 2;
    nq = size(M, 1);

    if nq < nr
        error('T:InvalidModel', ...
            'M must contain at least the two actuated rigid coordinates.');
    end

    validateattributes(M, {'numeric'}, ...
        {'2d','square','size',[nq nq],'finite'}, mfilename, 'M');
    validateattributes(D, {'numeric'}, ...
        {'2d','size',[nq nq],'finite'}, mfilename, 'D');
    validateattributes(K, {'numeric'}, ...
        {'2d','size',[nq nq],'finite'}, mfilename, 'K');
    validateattributes(H, {'numeric'}, ...
        {'column','numel',nq,'finite'}, mfilename, 'H');
    validateattributes(Q, {'numeric'}, ...
        {'2d','size',[2*nr 2*nr],'finite'}, mfilename, 'Q');
    validateattributes(R, {'numeric'}, ...
        {'2d','size',[nr nr],'finite'}, mfilename, 'R');

    % The inverse-dynamics part makes the rigid-coordinate plant behave
    % approximately as two decoupled double integrators:
    %       q_ddot_r = v
    % LQR therefore returns a commanded rigid-joint acceleration.
    A_lqr = [zeros(nr), eye(nr); zeros(nr), zeros(nr)];
    B_lqr = [zeros(nr); eye(nr)];

    % Reuse the gain while Q and R remain unchanged.  This avoids solving
    % the Riccati equation at every ODE evaluation.
    persistent Klqr Q_cached R_cached
    if isempty(Klqr) || ~isequal(Q,Q_cached) || ~isequal(R,R_cached)
        [Klqr,~,poles] = lqr(A_lqr, B_lqr, Q, R);
        if any(real(poles) >= 0)
            error('T:UnstableLQR', ...
                'The chosen Q and R do not produce a stable LQR design.');
        end
        Q_cached = Q;
        R_cached = R;
    end

    if numel(z) ~= nq || numel(z_dot) ~= nq
        error('T:InvalidStateDimension', ...
            'z and z_dot must have the same length as the augmented matrices.');
    end

    z = z(:);
    z_dot = z_dot(:);

    q_r = z(1:nr);
    q_dot_r = z_dot(1:nr);
    q_d_r = q_d(1:nr);
    q_dot_d_r = q_dot_d(1:nr);
    q_ddot_d_r = q_ddot_d(1:nr);

    x_error = [q_r - q_d_r; q_dot_r - q_dot_d_r];
    acceleration_command = q_ddot_d_r - Klqr*x_error;

    % Condense the unactuated flexible accelerations out of the nonlinear
    % equations while retaining their instantaneous dynamic effect.
    F_total = D*z_dot + K*z + H;

    if nq == nr
        M_equivalent = M;
        F_equivalent = F_total;
    else
        Mrr = M(1:nr, 1:nr);
        Mrf = M(1:nr, nr+1:end);
        Mfr = M(nr+1:end, 1:nr);
        Mff = M(nr+1:end, nr+1:end);

        Fr = F_total(1:nr);
        Ff = F_total(nr+1:end);

        M_equivalent = Mrr - Mrf*(Mff\Mfr);
        F_equivalent = Fr - Mrf*(Mff\Ff);
    end

    tau = M_equivalent*acceleration_command + F_equivalent;
    Tau1 = tau(1);
    Tau2 = tau(2);
end