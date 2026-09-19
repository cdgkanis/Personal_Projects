classdef Robot < handle
    properties
        DIM   = 2;
        N     = 2;
        NDoFs = 6;

        m_i   = []; % [-]
        L     = []; % [m]
        r     = []; % [m]
        M     = [];

        Mh    = []; %[Kg]
        Mp    =  0; %[Kg]
        Jh    = []; %[Kgm^2]
        Jp    =  0; %[Kgm^2]

        % Material
        E     =  0; %[Pa]
        ro    =  0; %[Kg/m^3]
        %nu    =  0; %[-]

        End_effector = struct();
        Controler = struct();

        f_bar  = [];
        f_sd   = [];
        N_mes = 0;
    end

    properties(Access = protected)
        % Calculated Characteristics
        d         = []; %[m]   Half Length of every link
        ro_linear = []; % Linear Density of every link

        I         = []; %[m^4]   Surface Inertia of every link
        Jo_tip    = []; %[Kgm^2] Mass Inertia of every link, regarding TIP
        Jo        = []; %[Kgm^2] Mass Inertia of every link, regarding COM

        % Generalized Coordinates
        q     = [];
        q_dot = [];

        u     = [];
        z     = [];

        % End - Properties
        ML        = [];
        JL        = [];
        MD        = [];

        % Eigen Frequencies
        beta      = [];


        % Mode Coeficients
        C         = [];

        % Matrices
        K         = [];
        D         = [];
        B         = [];
        H         = [];

        omega = [];
    end

    methods
        function Robot = Robot(Robot_params, Material_params)
            if (Robot_params.DIM >=2)
                Robot.DIM = Robot_params.DIM;
            else
                error("Robot: DIM has to be 2 or 3");
            end

            if (Robot_params.m_i >= 0)
                Robot.m_i = Robot_params.m_i;
            else
                error("Robot: Mass of Links has to be non-negative");
            end

            if (Robot_params.L > 0)
                Robot.L   = Robot_params.L;
                Robot.d   = Robot.L/2;
            else
                error("Robot: Length of Links has to be positive");
            end

            if (Robot_params.r >= 0)
                Robot.r   = Robot_params.r;
            else
                error("Robot: Equivalent radious of Links has to be positive");
            end

            Robot.Mh  = Robot_params.Mh;
            Robot.Mp  = Robot_params.Mp;
            Robot.Jh  = Robot_params.Jh;
            Robot.Jp  = Robot_params.Jp;

            Robot.E = Material_params.E;
            Robot.ro = Material_params.ro;
            %Robot.nu = Material_params.nu;

            % Calculating Resulting Characteristics
            Robot.N     = length(Robot.L);
            Robot.NDoFs = Robot.N + (Robot.DIM - 1)*sum(Robot.m_i);

            Robot.M   = Material_params.ro *pi* Robot.r .^2 .* Robot.L;%Robot_params. M_i

            Robot.ro_linear = Robot.M ./ Robot.L;
            Robot.Jo = zeros(Robot.N, 1);
            Robot.Jo(:) = computeInetriaCOM(Robot.L(:), Robot.r(:), Robot.M(:));
            Robot.I = zeros(Robot.N, 1);
            Robot.I(:)  = computeSurfaceInertia(Robot.r(:));

            Robot.q     = zeros(Robot.N + sum(Robot.m_i),1);
            Robot.q_dot = zeros(Robot.N + sum(Robot.m_i),1);

            Robot.C    = zeros(2, Robot.N, max(Robot.m_i));
            Robot.beta = zeros(   Robot.N, max(Robot.m_i));
            Robot.f_bar = zeros(  Robot.NDoFs+2,1);
            Robot.f_sd = zeros(  Robot.NDoFs+2,1);

            Robot.End_effector = Robot_params.End_effector;
            Robot.Controler = Robot_params.Controler;
        end

        function Robot = setState(Robot, q, q_dot)
            if length(q) == length(q_dot)
                if length(q) == sum(Robot.m_i) + Robot.N
                    Robot.q = q;
                    Robot.q_dot = q_dot;
                end
            end
        end

        function Print_properties(Robot)

            fprintf('\n========== SYSTEM PROPERTIES ==========\n');

            % General
            fprintf('\n-- General --\n');
            fprintf('Number of links (N): %d\n', Robot.N);

            fprintf('Number of Assumed Frequencies of links (m_i): [');
            fprintf('%d ', Robot.m_i);
            fprintf(']\n');

            % Dimensions
            fprintf('\n-- Dimensions [m] --\n');

            fprintf('Lengths (L):          [');
            fprintf('%.3f ', Robot.L);
            fprintf('] [m]\n');

            fprintf('COM Position (d):     [');
            fprintf('%.3f ', Robot.d);
            fprintf('] [m]\n');

            % Mass
            fprintf('\n-- Mass [kg] --\n');

            fprintf('Link mass (M):        [');
            fprintf('%.3f ', Robot.M);
            fprintf('] [kg]\n');

            fprintf('Hub mass (M_h):       [');
            fprintf('%.3f ', Robot.Mh);
            fprintf('] [kg]\n');

            fprintf('Payload mass (M_p): %.3f kg\n', Robot.Mp);

            % Inertia
            fprintf('\n-- Inertia [kg·m^2] --\n');

            fprintf('Link inertia (Jo):    [');
            fprintf('%.3e ', Robot.Jo);
            fprintf('] [kg·m^2]\n');

            fprintf('Hub inertia (Jh):     [');
            fprintf('%.3e ', Robot.Jh);
            fprintf('] [kg·m^2]\n');

            fprintf('Payload inertia (Jp): %.3e kg·m^2\n', Robot.Jp);

            % Structural Properties
            EI = Robot.E .* Robot.I;

            fprintf('\n-- Structural Properties --\n');

            fprintf('Density (ro): %.3e kg/m^3\n', Robot.ro);

            fprintf('Bending stiffness (EI): [');
            fprintf('%.3e ', EI);
            fprintf('] [N·m^2]\n');

            fprintf('\n=======================================\n\n');

        end

        function Print_End_properties(Robot)
            fprintf('\n======= Link End Characteristics =======\n');

            n = Robot.N;

            % --- ML ---
            fprintf('\n-- ML [kg] --\n');
            if ~isempty(Robot.ML)
                for i = 1:n
                    fprintf('Link %d: ML[%d] = %.4f\n', i, i, Robot.ML(i));
                end
            else
                fprintf('Print_End_properties: ML not computed yet.\n');
            end

            % --- JL ---
            fprintf('\n-- JL [kg·m^2] --\n');
            if ~isempty(Robot.JL)
                for i = 1:n
                    fprintf('Link %d: JL[%d] = %.6f\n', i, i, Robot.JL(i));
                end
            else
                fprintf('Print_End_properties: JL not computed yet.\n');
            end

            % --- MD ---
            fprintf('\n-- MD [kgm] --\n');
            if ~isempty(Robot.MD)
                for i = 1:n
                    fprintf('Link %d: ', i);
                    fprintf('%.4f ', Robot.MD(i));
                    fprintf('\n');
                end
            else
                fprintf('Print_End_properties: MD not computed yet.\n');
            end

            fprintf('========================================\n\n');
        end

        function Print_frequency_properties(obj)
            fprintf('\n======= Link Frequency Characteristics =======\n');

            n = obj.N;
            noAM = obj.m_i;
            % --- Beta ---
            fprintf('\n-- β [1/m] --\n');
            if ~isempty(obj.beta)
                for i = 1:n
                    for j = 1:noAM(i)
                        fprintf('Link %d: beta[%d, %d] = %.4f\n', i, i, j, obj.beta(i, j));
                    end
                end
            else
                fprintf('β not computed yet.\n');
            end

            % --- Omega ---
            fprintf('\n-- ω [rad/s] --\n');
            if ~isempty(obj.omega)
                for i = 1:n
                    for j = 1:noAM(i)
                        fprintf('Link %d: ω[%d, %d] = %.6f\n', i, i, j, obj.omega(i, j));
                    end
                end
            else
                fprintf('ω not computed yet.\n');
            end
            fprintf('========================================\n\n');

            % --- Frequency ---
            fprintf('\n-- f [Hz] --\n');
            if ~isempty(obj.omega)
                for i = 1:n
                    for j = 1:noAM(i)
                        fprintf('Link %d: f[%d, %d] = %.6f\n', i, i, j, obj.omega(i, j)/2/pi);
                    end
                end
            else
                fprintf('f not computed yet.\n');
            end
            fprintf('========================================\n\n');
        end

        % Give φ_ij values
        function phi = Phi(obj, i, j, x)
            C1 = obj.C(1,i,j);
            C2 = obj.C(2,i,j);

            beta_ij = obj.beta(i,j);
            phi = C1*sin(beta_ij*x) + C2*cos(beta_ij*x) - C1*sinh(beta_ij*x) - C2*cosh(beta_ij*x);
        end

        function phi = Phi_t(obj, i, j, x)
            C1 = obj.C(1,i,j);
            C2 = obj.C(2,i,j);

            beta_ij = obj.beta(i,j);
            phi = beta_ij*(C1*cos(beta_ij*x) - C1*cosh(beta_ij*x) - C2*sin(beta_ij*x) - C2*sinh(beta_ij*x));
        end

        function phi = Phi_tt(obj, i, j, x)
            C1 = obj.C(1,i,j);
            C2 = obj.C(2,i,j);

            beta_ij = obj.beta(i,j);
            phi = beta_ij^2*(- C1*sin(beta_ij*x) - C2*cos(beta_ij*x) - C1*sinh(beta_ij*x) - C2*cosh(beta_ij*x));
        end

        function Plot_Phi(obj, i, j)
            % Plot the mode shape Phi for link i, mode j

            Len = obj.L(i);                  % Link length
            x = linspace(0, Len, 1000);      % Spatial points along the link

            % Evaluate the mode shape
            phi = obj.Phi(i, j, x);

            % Create figure
            figure('Color', 'w');          % White background for publication

            % Plot mode shape with nice color
            plot(x, phi, 'Color', [0.0, 0.6, 0.6], 'LineWidth', 2);
            hold on;

            % Enhance appearance
            grid on;
            box on;

            ax = gca;
            ax.FontSize = 12;
            ax.LineWidth = 1;

            % Labels with clear formatting
            xlabel('Link Position x [m]', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel(['Mode Shape \Phi_{', num2str(i), ',', num2str(j), '} [-]'], 'FontSize', 12, 'FontWeight', 'bold');

            % Title
            title(['Mode Shape \Phi_{', num2str(i), ',', num2str(j), '}'], ...
                'FontSize', 14, 'FontWeight', 'bold');

            % Optional: add line at x=0 for reference
            plot([0, Len], [0, 0], '--k', 'LineWidth', 1);  % baseline at phi=0

            % --- Display C1 and C2 in a small textbox ---
            C1 = obj.C(1, i, j);
            C2 = obj.C(2, i, j);

            txt = sprintf('C_1 = %.4f\nC_2 = %.4f', C1, C2);

            % Place the textbox at 70% of x and 80% of y axis range
            xPos = 0.8*Len;
            yPos = 1.2*(max(phi)-min(phi)) + min(phi);

            text(xPos, yPos, txt, 'BackgroundColor', 'w', ...
                'EdgeColor', 'k', 'FontSize', 11, 'FontWeight', 'bold', ...
                'Margin', 5);


            hold off;
        end

        function Plot_Phi_link(obj, i)
            % Plot the mode shape Phi for link i, mode j

            Len = obj.L(i);                  % Link length
            x = linspace(0, Len, 1000);      % Spatial points along the link

            % Evaluate the mode shape
            phi1 = obj.Phi(i, 1, x);
            phi2 = obj.Phi(i, 2, x);

            % L_1_test = trapz(phi1.^2,x)
            % L_2_test = trapz(phi2.^2,x)


            % Create figure
            figure('Color', 'w');          % White background for publication

            % Plot mode shape with nice color
            plot(x, phi1, 'Color', [1, 0.6, 0.6], 'LineWidth', 2);
            hold on;
            plot(x, phi2, 'Color', [0, 1, 0.8], 'LineWidth', 2);
            hold on;
            legend(['Mode Shape \Phi_{', num2str(i),',',  num2str(1),'}'],...
                ['Mode Shape \Phi_{', num2str(i),',',  num2str(2),'}']);
            % Enhance appearance
            grid on;
            box on;

            ax = gca;
            ax.FontSize = 12;
            ax.LineWidth = 1;

            % Labels with clear formatting
            xlabel('Link Position x [m]', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel(['Mode Shape \Phi_{', num2str(i), ',', num2str(1), '}|\Phi_{', num2str(i), ',', num2str(2), '} [-]'], 'FontSize', 12, 'FontWeight', 'bold');

            % Title
            title(['Mode Shape \Phi_{', num2str(i), ',', num2str(1), '}|\Phi_{', num2str(i), ',', num2str(2), '}'], ...
                'FontSize', 14, 'FontWeight', 'bold');

            % Optional: add line at x=0 for reference
            plot([0, Len], [0, 0], '--k', 'LineWidth', 1);  % baseline at phi=0
            hold off;
        end

        % Give y_i
        function result = y(obj, i, x)
            noAM = obj.m_i(i);
            result = 0;

            for j= 1:noAM
                k = obj.N + sum(noAM(1:i-1)) + j;
                dij = obj.q(k);
                result = result + obj.Phi(i, j, x)*dij;
            end
        end

        function result = y_t(obj, i, x)
            noAM = obj.m_i(i);
            result = 0;

            for j= 1:noAM
                k = obj.N + sum(noAM(1:i-1)) + j;
                dij = obj.q(k);
                result = result + obj.Phi_t(i, j, x)*dij;
            end
        end

        function result = y_tt(obj, i, x)
            noAM = obj.m_i(i);
            result = 0;

            for j= 1:noAM
                k = obj.N + sum(noAM(1:i-1)) + j;
                result = result + obj.Phi_tt(i, j, x);
            end
        end

        function result = y_dot(obj, i, x)
            noAM = obj.m_i(i);
            result = 0;

            for j= 1:noAM
                k = obj.N + sum(noAM(1:i-1)) + j;
                dij_dot = obj.q_dot(k);
                result = result + obj.Phi(i, j, x)*dij_dot;
            end
        end

        function result = y_ddot(obj, i, x, q_ddot)
            noAM = obj.m_i(i);
            result = 0;

            for j= 1:noAM
                k = obj.N + sum(noAM(1:i-1)) + j;
                dij_ddot = q_ddot(k);
                result = result + obj.Phi(i, j, x)*dij_ddot;
            end
        end

        function result = y_t_dot(obj, i, x)
            noAM = obj.m_i(i);
            result = 0;

            for j= 1:noAM
                k = obj.N + sum(noAM(1:i-1)) + j;
                dij_dot = obj.q_dot(k);
                result = result + obj.Phi_t(i, j, x)*dij_dot;
            end
        end

        function result = y_t_ddot(obj, i, x, q_ddot)
            noAM = obj.m_i(i);
            result = 0;

            for j= 1:noAM
                k = obj.N + sum(noAM(1:i-1)) + j;
                dij_dot = q_ddot(k);
                result = result + obj.Phi_t(i, j, x)*dij_dot;
            end
        end

        function result = ys(obj, i, x, q)
            noAM = obj.m_i;
            result = 0;

            for j= 1:noAM(i)
                k = obj.N + sum(noAM(1:i-1)) + j;
                dij = q(k);
                result = result + obj.Phi(i, j, x)*dij;
            end
        end

        function result = ys_t(obj, i, x, q)
            noAM = obj.m_i;
            result = 0;

            for j= 1:noAM(i)
                k = obj.N + sum(noAM(1:i-1)) + j;
                dij = q(k);
                result = result + obj.Phi_t(i, j, x)*dij;
            end
        end

        function result = ys_dot(obj, i, x, q_dot)
            noAM = obj.m_i;
            result = 0;

            for j= 1:noAM(i)
                k = obj.N + sum(noAM(1:i-1)) + j;
                dij_dot = q_dot(k);
                result = result + obj.Phi(i, j, x)*dij_dot;
            end
        end

        function result = ys_t_dot(obj, i, x, q_dot)
            noAM = obj.m_i;
            result = 0;

            for j= 1:noAM(i)
                k = obj.N + sum(noAM(1:i-1)) + j;
                dij_dot = q_dot(k);
                result = result + obj.Phi_t(i, j, x)*dij_dot;
            end
        end

        function [y, y_dot, A, E, W_hat_] = FW_Kinematics_2D(obj, q, q_dot, link, point)
            S = [0 -1;1 0];

            A = zeros(2,2,2);
            E = zeros(2,2,2);
            W_hat_ = zeros(2,2);

            W_hat = eye(2);
            W_hat_dot = zeros(2, 2);

            ri = zeros(2,1);
            ri_dot = zeros(2,1);

            theta_i = q(1:obj.N);
            omega_i = q_dot(1:obj.N);

            for i = 1:obj.N

                if i > link
                    break;
                end

                A_i = [cos(theta_i(i)) -sin(theta_i(i));sin(theta_i(i)) cos(theta_i(i))];
                A_i_dot = S*A_i*omega_i(i);

                A(i,:,:) = A_i(:,:);

                W = W_hat*A_i;
                W_dot = W_hat_dot*A_i + W_hat*A_i_dot;

                i_pi = [point obj.ys(link,point,q)]';
                i_pi_dot = [0 obj.ys_dot(i,point,q_dot)]';

                y = ri + W*i_pi;
                y_dot = ri_dot + W_dot*i_pi + W*i_pi_dot;

                i_ri1     = [obj.L(i) obj.ys(i,obj.L(i),q)]';
                i_ri1_dot = [0 obj.ys_dot(i,obj.L(i),q_dot)]';

                ri     = ri     + W*i_ri1;
                ri_dot = ri_dot + W_dot*i_ri1 + W*i_ri1_dot;

                E_i = [1 -obj.ys_t(i,obj.L(i),q);obj.ys_t(i,obj.L(i),q) 1];
                E_i_dot = S*obj.ys_t_dot(i,obj.L(i),q_dot);

                E(i,:,:) = E_i(:,:);

                W_hat = W*E_i;
                W_hat_(:,:) = W_hat(:, :);
                W_hat_dot = W_dot*E_i + W*E_i_dot;
            end
        end

        function D = getDistanceMatrixHUB(obj)
            % GETDISTANCEMATRIX Computes an (N+1)x(N+1) matrix containing the
            % Euclidean distances between all nodes (Base, Joints, and Tip),
            % accounting for both rigid rotation and elastic deflection.

            nNodes = obj.N + 1;
            positions = zeros(nNodes, obj.DIM);

            % Zero vectors for velocity arguments required by your FW_Kinematics
            q_dot_zero = zeros(size(obj.q_dot));

            % 1. Position of Node 1 (Base is at the origin)
            % (Or call FW_Kinematics(obj.q, q_dot_zero, 1, 0) if base shifts)
            positions(:, :) = zeros(obj.N+1, obj.DIM);

            % 2. Collect absolute positions of all subsequent nodes
            for i = 1:obj.N
                % Node i+1 is located at the end (point = obj.L(i)) of link i
                [y_pos, ~] = obj.FW_Kinematics_2D(obj.q, q_dot_zero, i, obj.L(i));

                % Store only the relevant coordinates (X, Y up to obj.DIM)
                positions(i+1, :) = y_pos(1:obj.DIM)';
            end

            % 3. Compute the pairwise Euclidean distance matrix
            % pdist2 handles the matrix generation instantly without loops
            D = pdist2(positions, positions);
        end

        function Dcom = getDistanceMatrixCOM(obj)
            % GETDISTANCEMATRIXCOM Computes the Euclidean distance between every node
            % (Base + Joints) and the COM of every link.

            nNodes = obj.N + 1;
            nLinks = obj.N;

            % Positions of Base + Joints
            nodePos = zeros(nNodes,obj.DIM);

            % Positions of the COM of each link
            comPos = zeros(nLinks,obj.DIM);

            q_dot_zero = zeros(size(obj.q_dot));

            %-------------------------------------------------------------
            % Base
            %-------------------------------------------------------------
            nodePos(1,:) = zeros(1,obj.DIM);

            %-------------------------------------------------------------
            % Joint positions
            %-------------------------------------------------------------
            for i = 1:obj.N

                [y,~] = obj.FW_Kinematics_2D(obj.q,q_dot_zero,i,obj.L(i));

                nodePos(i+1,:) = y(1:obj.DIM)';

            end

            %-------------------------------------------------------------
            % COM positions
            %-------------------------------------------------------------
            for i = 1:obj.N

                [y,~] = obj.FW_Kinematics_2D(obj.q,q_dot_zero,i,obj.d(i));

                comPos(i,:) = y(1:obj.DIM)';

            end

            %-------------------------------------------------------------
            % Distance matrix
            % Row    -> node
            % Column -> COM of link
            %-------------------------------------------------------------
            Dcom = pdist2(nodePos,comPos);

        end

        function obj = Calculate_u(obj, i, j)
            ro_c = obj.ro*pi*obj.r(i)^2;
            L_c = obj.L(i);

            x = linspace(0, L_c, 10000);
            phi = obj.Phi(i, j, x);
            result = trapz( x, ro_c*phi);
            obj.u(i ,j) = result;
        end

        function obj = Calculate_z(obj, i, j, k)
            ro_c = obj.ro*pi*obj.r(i)^2;
            L_c = obj.L(i);

            x = linspace(0, L_c, 1000);
            phi_ij = obj.Phi(i, j, x);
            phi_ik = obj.Phi(i, k, x);

            obj.z(i, j, k) = trapz(x, ro_c*phi_ij.*phi_ik);
        end

        % Calculate ML_i | JL_i | MD_i
        function obj = Calculate_End_Characteristics(obj)

            n = obj.N;
            theta = obj.q(1:n);
            delta = obj.q(n+1:end);

            % ML
            obj.ML = zeros(n, 1);
            obj.ML(n) = obj.Mp;
            for i = n-1:-1:1
                % Sum link masses including payload
                obj.ML(i) = obj.ML(i+1) + obj.M(i+1);

                % Sum hub masses
                obj.ML(i) = obj.ML(i) + obj.Mh(i+1);
            end

            % --- MD (First Moment of Inertia (Unbalance) ) ---

            obj.MD = zeros(n,1);

            % % Link n
            % obj.MD(n)=0;
            %
            % % Link n-1
            % MD_i = obj.M(n)*obj.d(n)*cos(theta(n));
            % MD_i = MD_i + obj.Mp*obj.L(n)*cos(theta(n));
            %
            % % Compliance contribution
            % for j = 1:obj.m_i(n)
            %     obj.Calculate_u(n, j);
            %     u_ij = obj.u(n, j);
            %     idx = sum(obj.m_i(1:n-1)) + j;
            %
            %     MD_i = MD_i - u_ij*delta(idx)*sin(theta(n));
            % end
            %
            % % Hub contribution
            % for j = 1:obj.m_i(n)
            %     Phi_ij_e = obj.Phi(n, j, obj.L(n));
            %     idx = sum(obj.m_i(1:n-1)) + j;
            %
            %     MD_i = MD_i - obj.Mp*Phi_ij_e*delta(idx)*sin(theta(n));
            % end
            %
            % % Load
            %
            % obj.MD(n - 1) = MD_i;
            %
            % for i = n-2:-1:1
            %     % Rigid contribution
            %     MD_i = obj.MD(i+1) + obj.M(i+1)*obj.d(i+1)*cos(theta(i+1));
            %
            %     % Compliance contribution
            %     for j = 1:obj.m_i(i+1)
            %         obj.Calculate_u(i+1, j);
            %         u_ij = obj.u(i+1, j);
            %         idx = sum(obj.m_i(1:i)) + j;
            %
            %         MD_i = MD_i - u_ij*delta(idx)*sin(theta(i+1));
            %     end
            %
            %     % Hub contribution
            %     for j = 1:obj.m_i(i+1)
            %         Phi_ij_e = obj.Phi(i+1, j, obj.L(i+1));
            %         idx = sum(obj.m_i(1:i)) + j;
            %
            %         MD_i = MD_i - obj.Mh(i+2)*Phi_ij_e*delta(idx)*sin(theta(i+1));
            %     end
            %
            %     obj.MD(i) = MD_i;
            % end

            % --- JL (Second Moment of Inertia ) ---

            obj.JL = zeros(n,1);
            D_COM = obj.getDistanceMatrixCOM();
            D_HUB = obj.getDistanceMatrixHUB();

            % Link n
            obj.JL(n) = obj.Jp;

            % Link n-1
            JL_i =  obj.Jo(n) + obj.M(n)*D_COM(n,n)^2;
            % Hub
            JL_i = JL_i + obj.Jh(n);
            % Load
            JL_i = JL_i + obj.Jp + obj.Mp*D_HUB(n,n+1)^2;
            obj.JL(n-1) = JL_i;


            for i = n-2:-1:1
                JL_i = 0;

                for k = i+1:n
                    % Link
                    JL_i = JL_i + obj.Jo(k) +  obj.M(k)*D_COM(i+1,k)^2;

                    % Hub
                    JL_i = JL_i + obj.Jh(k) + obj.Mh(k)*D_HUB(i+1,k)^2;
                end

                % Load
                JL_i = JL_i + obj.Jp + obj.Mp*D_HUB(i+1,n)^2;
                obj.JL(i) = JL_i;
            end
        end

        function F = Construct_F_2D(obj, i, beta_ij)
            L_ = obj.L(i);
            b = beta_ij;

            EI_c = obj.E* obj.I(i);

            omega_ij = beta_ij^2*sqrt( EI_c/ obj.ro/(pi*obj.r(i)^2));

            ML_c = obj.ML(i);
            JL_c = obj.JL(i);
            MD_c = obj.MD(i);

            der = -1;

            % --- Coefficients from boundary conditions ---
            F11 =     EI_c           *b^2*(-sin(b*L_) - sinh(b*L_)) ...
                + JL_c*omega_ij^2*b  *( cos(b*L_) - cosh(b*L_))*(der)...
                + MD_c*omega_ij^2*    ( sin(b*L_) - sinh(b*L_))*(der);

            F12 =     EI_c           *b^2*(-cos(b*L_) - cosh(b*L_)) ...
                + JL_c*omega_ij^2*b  *(-sin(b*L_) - sinh(b*L_))*(der)...
                + MD_c*omega_ij^2    *( cos(b*L_) - cosh(b*L_))*(der);

            F21 =     EI_c           *b^3*(-cos(b*L_) - cosh(b*L_)) ...
                - ML_c*omega_ij^2*    ( sin(b*L_) - sinh(b*L_))*(der)...
                - MD_c*omega_ij^2*b  *( cos(b*L_) - cosh(b*L_))*(der);

            F22 =     EI_c           *b^3*( sin(b*L_) - sinh(b*L_)) ...
                - ML_c*omega_ij^2    *( cos(b*L_) - cosh(b*L_))*(der)...
                - MD_c*omega_ij^2*b  *(-sin(b*L_) - sinh(b*L_))*(der);

            F = [F11 F12;F21 F22];
        end

        function obj = Calculate_Eigen_Frequencies(obj, i)

            nModes = obj.m_i(i);

            solsetings = optimset('TolX', 1e-10);

            char_eq = @(beta) min(svd(obj.Construct_F_2D(i,beta)));

            %------------------------------------------------------
            % Scan
            %------------------------------------------------------
            betaMax = 80/obj.L(i);
            Nscan   = 5000;

            beta_ = linspace(0.01,betaMax,Nscan);

            f = zeros(size(beta_));

            for k = 1:Nscan
                f(k) = char_eq(beta_(k));
            end

            %------------------------------------------------------
            % Find candidate minima
            %------------------------------------------------------
            idx = islocalmin(f);
            idx = find(idx);

            roots = [];

            %------------------------------------------------------
            % Refine using direct minimization (NOT fzero)
            %------------------------------------------------------
            for k = 1:length(idx)

                if idx(k) == 1 || idx(k) == Nscan
                    continue;
                end

                b1 = beta_(idx(k)-1);
                b2 = beta_(idx(k)+1);

                % local minimization instead of root finding
                r_ = fminbnd(char_eq, b1, b2, solsetings);

                % avoid duplicates
                if isempty(roots) || all(abs(roots - r_) > 1e-4)
                    roots(end+1) = r_;
                end
            end

            roots = sort(roots);

            %------------------------------------------------------
            % Store
            %------------------------------------------------------
            nFound = min(nModes, length(roots));

            for j = 1:nFound
                obj.beta(i,j) = roots(j);

                obj.omega(i,j) = roots(j)^2 * ...
                    sqrt(obj.E*obj.I(i)/(obj.ro*pi*obj.r(i)^2));
            end
        end

        function obj = Calculate_C(obj, i, j)
            % --- Coefficients from boundary conditions ---
            betta = obj.beta(i, j);
            L_ = obj.L(i);

            F = obj.Construct_F_2D(i, betta);

            % Choose the more stable equation
            if abs(F(1,1)) < 1e-14 && abs(F(2,1)) < 1e-14
                error('Singular boundary system');
            end
            if abs(F(1, 2)) > abs(F(2, 2))
                Bi = -F(1, 2) / F(1, 1);
            else
                Bi = -F(2, 2) / F(2, 1);
            end

            if abs(betta) <= 1e-14
                error('beta too small');
            end
            if abs(Bi) < 1e-12
                warning('Bi ~ 0 : mode dominated by second basis function');
            end

            C1 = Bi;
            C2 = 1;

            x = linspace(0, L_, 10000);
            %Phi
            Phi_ij= (- C1*sin(betta*x) - C2*cos(betta*x) - C1*sinh(betta*x) - C2*cosh(betta*x));

            zijj = trapz(x,Phi_ij.^2);
            norm = -sqrt(L_/abs(zijj));

            C1 = C1*norm;
            C2 = C2*norm;
            % --- Store coefficients ---

            obj.C(1,i,j) = C1;
            obj.C(2,i,j) = C2;
        end

        function obj = Find_C_Coeficients(obj, i)
            noam_i = obj.m_i(i);
            for j = 1:noam_i
                obj.Calculate_C(i, j);
            end
        end

        function Check(obj)
            for i = 1:obj.N
                ml = obj.ML(i);
                md = obj.MD(i);
                jl = obj.JL(i);
                L_ = obj.L(i);

                Orth = zeros(obj.m_i(i), obj.m_i(i));

                for j = 1:obj.m_i(i)
                    Phii1e   = obj.Phi  (i, j, L_);
                    Phii1e_t = obj.Phi_t(i, j, L_);
                    for k = j+1:obj.m_i(i)
                        Phii2e   = obj.Phi  (i, k, L_);
                        Phii2e_t = obj.Phi_t(i, k, L_);

                        obj.Calculate_z(i,j,k);
                        Orth(j, k) = [Phii1e Phii1e_t]*[ml 0.5*md;0.5*md jl]*[Phii2e;Phii2e_t] + obj.z(i,j,k);
                    end
                end
                disp("Orthogonality");
                disp(Orth);
            end
        end

        % Matrices
        function Calculate_K_matrix(obj)
            NoL    = obj.N;
            NoDoFs = obj.NDoFs;
            NoAM   = obj.m_i;

            obj.K = zeros(NoDoFs, NoDoFs);

            % for i = 1:NoL
            %     for j = 1:NoAM(i)
            %         idx = NoL + sum(NoAM(1:i-1)) + j;
            %
            %         obj.K(idx, idx) = obj.M(i)*obj.omega(i,j)^2;
            %     end
            % end


            for i = 1:NoL
                L_c = obj.L(i);
                x = linspace(0, L_c, 10000);
                y_tt   = obj.y_tt(i,x);

                for j = 1:NoAM(i)
                    k1 = obj.N + sum(NoAM(1:i-1)) + j;
                    phi_ttj = obj.Phi_tt(i, j, x);
                    for g= 1:NoAM(i)
                        k2 = obj.N + sum(NoAM(1:i-1)) + g;
                        phi_ttg = obj.Phi_tt(i, g, x);

                        result = trapz( x, phi_ttj.*phi_ttg);
                        obj.K(k1, k2) = obj.E*obj.I(i)*result;
                    end
                end
            end
        end

        function Print_K_matrix(obj)
            disp("Stighness Matrix");
            disp(obj.K);
        end

        function Calculate_D_matrix(obj)
            NoDoFs = obj.NDoFs;
            for i = 1:NoDoFs
                if obj.K(i,i) >0
                    obj.D(i, i) = 0.0001*sqrt(obj.K(i, i));
                else
                    obj.D(i, i) = 0.0001;
                end
            end
        end

        function Print_D_matrix(obj)
            disp("Dumping Matrix");
            disp(obj.D);
        end

        function Calculate_B_matrix(obj)

            nq = obj.NDoFs;
            obj.B  = zeros(nq,nq);

            [xi, wi] = gaussLegendre8();

            % ==========================================================
            % Distributed mass of every flexible link
            % ==========================================================
            for i = 1:obj.N

                Li   = obj.L(i);
                rhoi = obj.ro * pi * obj.r(i)^2;

                % Map Gauss points from [-1,1] to [0,Li]
                xg = 0.5 * Li * (xi + 1);
                wg = 0.5 * Li * wi;

                for k = 1:numel(xg)
                    Jv = obj.VelocityJacobian2D(i, xg(k), obj.q);
                    obj.B = obj.B + rhoi * wg(k) * (Jv.' * Jv);
                end
            end

            % ==========================================================
            % Hub masses
            % Mh(i) is located at the origin of link i. Therefore Mh(1)
            % is fixed at the base and, for i > 1, Mh(i) is located at
            % the end of link i-1.
            % ==========================================================
            for i = 2:obj.N

                if obj.Mh(i) ~= 0
                    Jv = obj.VelocityJacobian2D(i-1, obj.L(i-1), obj.q);
                    obj.B = obj.B + obj.Mh(i) * (Jv.' * Jv);
                end
            end

            % ==========================================================
            % Payload at tip of final link
            % ==========================================================
            if obj.Mp ~= 0
                Jv = obj.VelocityJacobian2D(obj.N, obj.L(obj.N), obj.q);
                obj.B  = obj.B + obj.Mp * (Jv.' * Jv);
            end

            % ==========================================================
            % Rotational inertia of hubs and payload
            % Under the small-deflection assumption, the absolute angle
            % of hub i contains all preceding link-end elastic slopes.
            % ==========================================================
            for i = 1:obj.N
                Jw = zeros(1,nq);
                Jw(1:i) = 1;

                for k = 1:i-1
                    for j = 1:obj.m_i(k)
                        idx = obj.N + sum(obj.m_i(1:k-1)) + j;
                        Jw(idx) = obj.Phi_t(k,j,obj.L(k));
                    end
                end

                obj.B = obj.B + obj.Jh(i) * (Jw.' * Jw);
            end

            Jwp = zeros(1,nq);
            Jwp(1:obj.N) = 1;

            for i = 1:obj.N
                for j = 1:obj.m_i(i)
                    idx = obj.N + sum(obj.m_i(1:i-1)) + j;
                    Jwp(idx) = obj.Phi_t(i,j,obj.L(i));
                end
            end

            obj.B = obj.B + obj.Jp * (Jwp.' * Jwp);

            % Numerical cleanup
            obj.B = 0.5*(obj.B + obj.B.');

        end

        function Jv = PositionJacobian2D(obj, link, point, epsq)

            nq = obj.NDoFs;
            q0 = obj.q;
            qdot0 = zeros(nq,1);

            Jv = zeros(2,nq);

            for k = 1:nq

                dq = zeros(nq,1);
                dq(k) = epsq;

                [rp, ~] = obj.FW_Kinematics_2D(q0 + dq, qdot0, link, point);
                [rm, ~] = obj.FW_Kinematics_2D(q0 - dq, qdot0, link, point);

                Jv(:,k) = (rp - rm)/(2*epsq);
            end
        end

        function Jv = VelocityJacobian2D(obj, link, point, q)

            nq = obj.NDoFs;
            Jv = zeros(2,nq);

            for k = 1:nq
                unitVelocity = zeros(nq,1);
                unitVelocity(k) = 1;

                [~, pointVelocity] = obj.FW_Kinematics_2D( ...
                    q, unitVelocity, link, point);

                Jv(:,k) = pointVelocity;
            end
        end

        function aBias = PointBiasAcceleration2D(obj, link, point, q, q_dot)

            if norm(q_dot,inf) == 0
                aBias = zeros(2,1);
                return;
            end

            % Directional derivative of J(q)*q_dot along q_dot.  The
            % perturbation is scaled so that no generalized coordinate is
            % displaced by more than approximately 1e-5.
            epsTime = 1e-5/max(1,norm(q_dot,inf));

            [~, velocityPlus] = obj.FW_Kinematics_2D( ...
                q + epsTime*q_dot, q_dot, link, point);
            [~, velocityMinus] = obj.FW_Kinematics_2D( ...
                q - epsTime*q_dot, q_dot, link, point);

            aBias = (velocityPlus - velocityMinus)/(2*epsTime);
        end

        function Calculate_H_matrix(obj)

            nq = obj.NDoFs;
            q0 = obj.q;
            qdot0 = obj.q_dot;

            if length(q0) ~= nq || length(qdot0) ~= nq
                error('Robot: state dimension is inconsistent with NDoFs');
            end

            % For a point mass with velocity v = J(q)*q_dot,
            %   B_point = m*J'*J
            %   H_point = m*J'*(J_dot*q_dot).
            % This is the same Lagrange/Christoffel vector as differentiating
            % the full B matrix, but it avoids 2*NDoFs complete mass-matrix
            % evaluations on every ODE call.
            obj.B = zeros(nq,nq);
            obj.H = zeros(nq,1);
            [xi, wi] = gaussLegendre8();

            % Distributed mass of every flexible link.
            for i = 1:obj.N
                Li = obj.L(i);
                rhoi = obj.ro*pi*obj.r(i)^2;
                xg = 0.5*Li*(xi + 1);
                wg = 0.5*Li*wi;

                for k = 1:numel(xg)
                    Jv = obj.VelocityJacobian2D(i,xg(k),q0);
                    aBias = obj.PointBiasAcceleration2D( ...
                        i,xg(k),q0,qdot0);

                    massWeight = rhoi*wg(k);
                    obj.B = obj.B + massWeight*(Jv.'*Jv);
                    obj.H = obj.H + massWeight*(Jv.'*aBias);
                end
            end

            % Translational kinetic energy of the moving hubs.
            for i = 2:obj.N
                if obj.Mh(i) ~= 0
                    Jv = obj.VelocityJacobian2D(i-1,obj.L(i-1),q0);
                    aBias = obj.PointBiasAcceleration2D( ...
                        i-1,obj.L(i-1),q0,qdot0);

                    obj.B = obj.B + obj.Mh(i)*(Jv.'*Jv);
                    obj.H = obj.H + obj.Mh(i)*(Jv.'*aBias);
                end
            end

            % Translational kinetic energy of the payload.
            if obj.Mp ~= 0
                Jv = obj.VelocityJacobian2D(obj.N,obj.L(obj.N),q0);
                aBias = obj.PointBiasAcceleration2D( ...
                    obj.N,obj.L(obj.N),q0,qdot0);

                obj.B = obj.B + obj.Mp*(Jv.'*Jv);
                obj.H = obj.H + obj.Mp*(Jv.'*aBias);
            end

            % Hub and payload angular-velocity Jacobians are constant for
            % the assumed-mode coordinates, so they contribute to B but
            % have zero rotational bias contribution to H.
            for i = 1:obj.N
                Jw = zeros(1,nq);
                Jw(1:i) = 1;

                for k = 1:i-1
                    for j = 1:obj.m_i(k)
                        idx = obj.N + sum(obj.m_i(1:k-1)) + j;
                        Jw(idx) = obj.Phi_t(k,j,obj.L(k));
                    end
                end

                obj.B = obj.B + obj.Jh(i)*(Jw.'*Jw);
            end

            Jwp = zeros(1,nq);
            Jwp(1:obj.N) = 1;

            for i = 1:obj.N
                for j = 1:obj.m_i(i)
                    idx = obj.N + sum(obj.m_i(1:i-1)) + j;
                    Jwp(idx) = obj.Phi_t(i,j,obj.L(i));
                end
            end

            obj.B = obj.B + obj.Jp*(Jwp.'*Jwp);
            obj.B = 0.5*(obj.B + obj.B.');

        end

        function Print_B_matrix(obj)
            disp("Mass Matrix");
            disp(obj.B);
        end

        function Print_H_matrix(obj)
            disp("Coriolis and Centrifugal Vector");
            disp(obj.H);
        end

        function [B, D, K, H] = getDynamicMatrices(obj)
            B = obj.B;
            D = obj.D;
            K = obj.K;
            H = obj.H;
        end

        function Mod_vec = applied_force_mod(obj, F_E)
            % Generalized force produced by an end-effector force.
            % Specialized for a two-link planar robot.
            %
            % F_E must be a 2x1 vector expressed in the end-effector frame.

            N_dofs = obj.NDoFs;
            Mod_vec = zeros(N_dofs,1);

            if obj.N ~= 2
                error('applied_force_mod is implemented only for two links.');
            end

            F_E = F_E(:);

            if numel(F_E) ~= 2
                error('F_E must be a two-component force vector.');
            end

            S = [0 -1;
                1  0];

            % Obtain individual link rotations and elastic rotations
            [~, ~, A, E_, ~] = obj.FW_Kinematics_2D( ...
                obj.q, obj.q_dot, 2, obj.L(2));

            % A and E_ are 2x2x2 arrays. Extract their matrix slices.
            A1 = squeeze(A(1,:,:));
            A2 = squeeze(A(2,:,:));

            E1 = squeeze(E_(1,:,:));

            % Frame matrices appearing in the derived equations
            W1 = A1;
            W_hat_1 = W1*E1;
            W2 = W_hat_1*A2;

            W = cell(2,1);
            W{1} = W1;
            W{2} = W2;

            l = obj.L;

            % Flexible transverse displacements at the link ends
            y1 = obj.ys(1, l(1), obj.q);
            y2 = obj.ys(2, l(2), obj.q);

            r1 = [l(1); y1];
            r2 = [l(2); y2];

            % Generalized force associated with theta_1
            Mod_vec(1) = F_E.'*W2.'*S*A1*(r1 + E1*r2);

            % Generalized force associated with theta_2
            Mod_vec(2) = F_E.'*W2.'*W_hat_1*S*A2*r2;

            % Generalized forces associated with delta_ij
            k = obj.N + 1;

            for i = 1:obj.N
                for j = 1:obj.m_i(i)

                    phi_ij = obj.Phi(i, j, l(i));

                    Mod_vec(k) = ...
                        F_E.'*W2.'*W{i}*[0; phi_ij];

                    if i == 1
                        phi_t_1j = obj.Phi_t(1, j, l(1));

                        Mod_vec(k) = Mod_vec(k) ...
                            + F_E.'*W2.'*W1*S*r2*phi_t_1j;
                    end

                    k = k + 1;
                end
            end
        end

        function r_n1_ddot = EEF_Base_acceleration(obj,q_ddot)
            l1 = obj.L(1);
            l2 = obj.L(2);

            y1_l1 = obj.y(1, l1);
            y2_l2 = obj.y(2, l2);
            y1_l1_dot = obj.y_dot(1, l1);
            y2_l2_dot = obj.y_dot(2, l2);
            y1_l1_ddot = obj.y_ddot(1, l1, q_ddot);
            y2_l2_ddot = obj.y_ddot(2, l2, q_ddot);

            [~, ~, A, E_] = obj.FW_Kinematics_2D(obj.q, obj.q_dot, 2, l2);

            S = [0 -1;
                1  0];

            % A and E_ are 2x2x2 arrays. Extract their matrix slices.
            A1 = squeeze(A(1,:,:));
            A1_dot  = S*A1*obj.q_dot(1);
            A1_ddot = S*(A1_dot*obj.q_dot(1) + A1*q_ddot(1));

            A2 = squeeze(A(2,:,:));
            A2_dot = S*A2*obj.q_dot(2);
            A2_ddot = S*(A2_dot*obj.q_dot(2) + A2*q_ddot(2));

            E1 = squeeze(E_(1,:,:));
            E1_dot = S*obj.y_t_dot(1,l1);
            E1_ddot = S*obj.y_t_ddot(1,l1, q_ddot);

            W1 = A1;
            W_hat_1 = W1*E1;
            W2 = W_hat_1*A2;

            W1_dot = A1_dot;
            W2_dot = (W1_dot*E1 + W1*E1_dot)*A2 + W_hat_1*A2_dot;

            W1_ddot = A1_ddot;
            W2_ddot = (W1_ddot*E1 + 2*W1_dot*E1_dot + W1*E1_ddot)*A2+ 2*(W1_dot*E1 + W1*E1_dot)*A2_dot + W_hat_1*A2_ddot;

            r_n1_ddot = W1_ddot*[l1;y1_l1] + 2*W1_dot*[0;y1_l1_dot] + W1*[0;y1_l1_ddot] +...
                + W2_ddot*[l2;y2_l2] + 2*W2_dot*[0;y2_l2_dot] + W2*[0;y2_l2_ddot] ;
        end

        function VEC2N = shape_functions(obj)

            VEC2N = zeros(2,obj.NDoFs);

            l(:) = obj.L(:);
            l1 = l(1);
            l2 = l(2);

            y1_l1 = obj.y(1, l1);
            y2_l2 = obj.y(2, l2);

            [~, ~, A, E_] = obj.FW_Kinematics_2D(obj.q, obj.q_dot, 2, l2);

            S = [0 -1;
                1  0];

            % A and E_ are 2x2x2 arrays. Extract their matrix slices.
            A1 = squeeze(A(1,:,:));
            dA1_dq  = S*A1;% d/d(theta_1) only

            A2 = squeeze(A(2,:,:));
            dA2_dq = S*A2;% d/d(theta_2) only

            E1 = squeeze(E_(1,:,:));

            W1 = A1;
            W_hat_1 = W1*E1;
            W2 = W_hat_1*A2;

            W = cell(2,1);
            W{1} = W1;
            W{2} = W2;

            VEC2N(:,1) = dA1_dq*[l1;y1_l1] + dA1_dq*E1*A2*[l2;y2_l2];% d/d(theta_1)
            VEC2N(:,2) = W_hat_1*dA2_dq*[l2;y2_l2];% d/d(theta_2)

            for i = 1:obj.N
                for j = 1:obj.m_i(i)
                    k = obj.N + sum(obj.m_i(1:i-1)) + j;

                    dy_ddij = obj.Phi(i, j, l(i));

                    VEC2N(:,k) = W{i}*[0;dy_ddij];

                    if i == 1
                        dE1_d1j = S*obj.Phi_t(1,j,l(1));

                        VEC2N(:,k) = VEC2N(:,k) + A1*dE1_d1j*A2*[l2;y2_l2];
                    end
                end
            end
        end
    
        function omg = renew_statistics(obj, omega_d)
            for i = 1:(obj.NDoFs+2)
            % renew statistics
                n = obj.N_mes;
                x = omega_d(i)/2/pi;

                mean_old = obj.f_bar(i);
                sd_old   = obj.f_sd(i);

                mean_new = (n*mean_old + x)/(n+1);

                if n == 0

                    sd_new = 0;

                elseif n == 1

                    % After second measurement
                    sd_new = sqrt( ...
                        (x - mean_old)*(x - mean_new) );

                else

                    var_new = ( ...
                        (n-1)*sd_old^2 + ...
                        (x - mean_old)*(x - mean_new) ) / n;

                    sd_new = sqrt(max(var_new,0));

                end

                obj.f_bar(i) = mean_new;
                obj.f_sd(i)  = sd_new;
            end
            obj.N_mes = obj.N_mes +1;
        end
    end
end

%% Helper Functions

function Jo = computeInetriaHUB(L, r, m)
Jo = 1/3*m.*L.*L + 1/4*m.*r.*r;
end

function Jo = computeInetriaCOM(L, r, m)
Jo =  1/12*m.*L.*L + 1/4*m.*r.*r;
end

function I = computeSurfaceInertia(r)
I = pi*r.^4/4;
end

function A = Zrotation(theta)
A = [cos(theta) -sin(theta) 0;...
    sin(theta) cos(theta) 0;...
    0 0 1];
end

function [xi,wi] = gaussLegendre8()
% Eight-point Gauss-Legendre rule on [-1,1].
xi = [-0.9602898564975363;...
    -0.7966664774136267;...
    -0.5255324099163290;...
    -0.1834346424956498;...
    0.1834346424956498;...
    0.5255324099163290;...
    0.7966664774136267;...
    0.9602898564975363];

wi = [0.1012285362903763;...
    0.2223810344533745;...
    0.3137066458778873;...
    0.3626837833783620;...
    0.3626837833783620;...
    0.3137066458778873;...
    0.2223810344533745;...
    0.1012285362903763];
end