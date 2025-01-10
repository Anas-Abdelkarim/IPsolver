function augmented_lagrangian_example
    % Example of implementing Algorithm 4.4 in MATLAB
    % Simple quadratic problem: minimize f(x) = 0.5*x'*Q*x + c'*x
    % subject to Ax <= b (convex constraint)

    % Problem data (adjust these as needed)
    Q = [3, 2; 2, 6];   % Hessian matrix
    c = [-2; -6];       % Linear term
    A = [1, 1];         % Inequality constraint
    b = 1;              % Right-hand side of constraint

    % Parameters
    rho = 1;            % Penalty parameter
    gamma = 1.5;        % Penalty scaling
    tau = 0.1;          % Feasibility tolerance
    max_iter = 100;     % Maximum iterations
    tol = 1e-6;         % Tolerance for stopping
    
    % Initialization
    x = [0; 0];         % Initial guess for x
    lambda = 0;         % Initial Lagrange multiplier
    B = [-10, 10];  % Boundaries for projection (box constraint)

    % Iterative optimization using Augmented Lagrangian
    for k = 1:max_iter
        % Step 1: Solve Augmented Lagrangian subproblem
        x = minimize_augmented_lagrangian(Q, c, A, b, x, lambda, rho);
        
        % Step 2: Update Lagrange multiplier
        Ax_plus_lambda_rho = A*x + lambda/rho;
        lambda_new = lambda + rho * (A*x - projection_box(Ax_plus_lambda_rho, B));
        
        % Step 3: Compute constraint violation and update rho
        V = norm(A*x - projection_box(Ax_plus_lambda_rho, B));
        if V > tol
            rho = gamma * rho;
        end
        
        % Step 4: Check for convergence
        if norm(A*x - b) < tol
            fprintf('Converged in %d iterations.\n', k);
            break;
        end
        
        lambda = lambda_new;  % Update multiplier
    end
    
    % Display result
    disp('Optimal solution:');
    disp(x);
end

% Function to minimize the augmented Lagrangian (quadratic in this case)
function x = minimize_augmented_lagrangian(Q, c, A, b, x0, lambda, rho)
    % Use MATLAB's quadprog to solve quadratic subproblem
    H = Q + rho * (A' * A);      % Augmented Hessian
    f = c - A' * (lambda - rho * b);  % Augmented linear term
    
    % Solve the quadratic subproblem
    options = optimoptions('quadprog', 'Display', 'none');
    x = quadprog(H, f, [], [], [], [], [], [], [], options);
end

% Function to project onto a box constraint
function z_proj = projection_box(z, B)
    % B is a matrix where each row is [lower_bound, upper_bound]
    z_proj = max(B(:, 1), min(z, B(:, 2)));
end
