close all
clearvars
clc

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Course of TOPOLOGY OPTIMISATION
%             A.Y. 2025/2026
%     Mechanical Engineering @ PoliMi
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -------------------------
% 1. LOAD DATA (Common for all cases)
% -------------------------
nodes = load('nodes.txt');
elements = load('elements.txt');
n_el = size(elements,1);

% Set Constraints (Fixed nodes at x=0)
constrNodes = find(nodes(:,1)==0);
constr = [sort(repmat(constrNodes, [2 1])) repmat([1;2], [length(constrNodes) 1])];

% Set External Load (6000N at node 49, -y direction)
ext_load = [49 2 -6000];

% Problem constants
x_lb = 1E-4*ones(n_el,1);           % Lower bound
x_ub = 800*ones(n_el,1);            % Upper bound
[V0_max, lengths] = computeVolume(x_ub, nodes, elements); % Max possible volume

% Initial Gradient of Constraint (Constant)
grad_constr_fixed = gradConstr(nodes, elements);

%% -------------------------
% 2. MAIN LOOP OVER ALPHA VALUES
% -------------------------
alpha_values = [0.05, 0.10, 0.15, 0.20];

for k = 1:length(alpha_values)
    alpha = alpha_values(k);
    
    % Print separator for clarity
    fprintf('\n========================================\n');
    fprintf('      RUNNING FOR ALPHA = %.2f\n', alpha);
    fprintf('========================================\n');
    
    % --- Initialize Problem for this Alpha ---
    V_target = alpha * V0_max;
    
    problem.x_ub = x_ub;
    problem.x_lb = x_lb;
    problem.V = V_target;
    problem.computeVolume = @computeVolume;
    problem.n_el = n_el;
    problem.nodes = nodes;
    problem.elements = elements;
    problem.grad_constr = grad_constr_fixed;
    problem.constr = constr;      % Needed for gradPhi if using fallback
    problem.ext_load = ext_load;  % Needed for gradPhi if using fallback

    % --- Initialization ---
    x_old = (V_target/sum(lengths)) * ones(n_el,1); % Initial guess
    lambda_old = 1;                                 % Initial lambda
    max_iter = 100;
    
    % Storage vectors
    of_vec = zeros(max_iter,1);
    V_vec = zeros(max_iter,1);
    x_vec = zeros(max_iter, n_el);
    lambda_vec = zeros(max_iter,1);

    % Initial Analysis
    K = buildStiffnessMatrix(nodes, elements, x_old);
    [u, F] = computeStructure(constr, ext_load, K);
    
    of_vec(1) = F'*u;
    V_vec(1) = computeVolume(x_old, nodes, elements);
    x_vec(1,:) = x_old';
    lambda_vec(1) = lambda_old;

    % --- CONLIN OPTIMIZATION LOOP ---
    for ii = 2:max_iter
        % 1. Analysis & Gradients
        K = buildStiffnessMatrix(nodes, elements, x_old);
        [u, F] = computeStructure(constr, ext_load, K);
        problem.u = u; % Store u for gradPhi to avoid re-calculation
        
        grad_of = gradOf(u, nodes, elements);
        problem.grad_of = grad_of;

        % 2. Solve Dual Problem (Find Lambda & x_new)
        % Using try-catch to prevent crash if fsolve fails
        try
            options = optimoptions('fsolve', 'Display', 'off');
            lambda_new = fsolve(@(lam) gradPhi(lam, x_old, problem), lambda_old, options);
        catch
            lambda_new = lambda_old; % Fallback
        end
        
        [~, x_new] = gradPhi(lambda_new, x_old, problem);

        % 3. Store Data
        lambda_vec(ii) = lambda_new;
        x_vec(ii,:) = x_new';
        
        % Re-evaluate for history storage
        K_new = buildStiffnessMatrix(nodes, elements, x_new);
        [u_new, F_new] = computeStructure(constr, ext_load, K_new);
        
        of_vec(ii) = F_new' * u_new;
        V_vec(ii) = computeVolume(x_new, nodes, elements);

        % 4. Update
        x_old = x_new;
        lambda_old = lambda_new;
    end
    
    x_opt = x_vec(max_iter,:)';

    % --- POST PROCESSING FOR THIS ALPHA ---
    
    % Compute final stresses
    K_opt = buildStiffnessMatrix(nodes, elements, x_opt);
    u_opt = computeStructure(constr, ext_load, K_opt);
    sigma_opt = computeStress(u_opt, x_opt, nodes, elements);
    
    % --- REPORT REQUESTED VALUES ---
    target_elements = [64; 77; 82];
    fprintf('Results for Alpha = %.2f:\n', alpha);
    fprintf('-------------------------------------------------\n');
    fprintf('| Elem | Area (mm^2)    | Stress (MPa)      |\n');
    fprintf('-------------------------------------------------\n');
    for i = 1:length(target_elements)
        el_idx = target_elements(i);
        fprintf('|  %2d  | %12.4f   | %15.4e   |\n', ...
            el_idx, x_opt(el_idx), sigma_opt(el_idx));
    end
    fprintf('-------------------------------------------------\n');

    % --- PLOTTING ---
    
    % Plot 1: Optimized Structure
    figure('Name', ['Structure alpha=' num2str(alpha)]);
    hold on
    for ii = 1:n_el
        % Calculate width for visualization
        w_vis = 0.5 + 5 * (x_opt(ii) - x_lb(ii)) / (max(x_opt) - x_lb(ii));
        
        plot([nodes(elements(ii,1), 1), nodes(elements(ii,2), 1)], ...
             [nodes(elements(ii,1), 2), nodes(elements(ii,2), 2)], ...
             'Color', 'k', 'LineWidth', w_vis);
    end
    % Draw nodes
    plot(nodes(:,1), nodes(:,2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
    axis equal; axis off;
    title(['Optimized Structure (\alpha = ' num2str(alpha) ')']);

    % Plot 2: Compliance History
    figure('Name', ['Compliance alpha=' num2str(alpha)]);
    plot(1:max_iter, of_vec, '-o', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
    grid on;
    xlabel('Iterations'); ylabel('Compliance [Nmm]');
    title(['Compliance History (\alpha = ' num2str(alpha) ')']);
    
    % Plot 3: Volume History
    figure('Name', ['Volume alpha=' num2str(alpha)]);
    plot(1:max_iter, V_vec, '-o', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
    hold on
    yline(V_target, '--r', 'Target Volume', 'LineWidth', 2);
    grid on;
    xlabel('Iterations'); ylabel('Volume [mm^3]');
    title(['Volume History (\alpha = ' num2str(alpha) ')']);
    
end

disp('All simulations completed.');
