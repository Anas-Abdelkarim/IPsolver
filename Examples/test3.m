% Parameters
rho = 1;            % Initial penalty parameter
lambda = 0;         % Initial multiplier for equality constraint
mu = [0, 0];        % Initial multipliers for inequality constraints
tol = 1e-6;         % Tolerance for stopping criteria
max_iter = 1000;     % Maximum number of iterations
x_k = [10; 10];     % Initial guess

for k = 1:max_iter
    % Define the objective function
    f0 = (x_k(1) - 50)^2 + (x_k(2) - 25)^2; % Cost function
    
    % Constraints
    h = x_k(1) + x_k(2) - 20;  % Equality constraint
    g1 = x_k(1) - 15;  % g1 <= 0
    g2 = 7 - x_k(2);   % g2 <= 0
    g = [g1; g2];
    
    % Augmented Lagrangian
    L_aug = f0 + lambda * h + (rho/2) * h^2 + ...
            sum(mu .* max(0, g)) + (rho/2) * sum(max(0, g).^2);
    
    % Gradient of the Augmented Lagrangian (symbolically computed)
    grad_Lx = 2 * (x_k(1) - 50) + lambda + rho * h + ...
              rho * max(0, g1) * (g1 > 0); % Partial derivative w.r.t x
    grad_Ly = 2 * (x_k(2) - 25) + lambda + rho * h + ...
              rho * max(0, g2) * (g2 > 0); % Partial derivative w.r.t y
    
    % Update step using gradient descent
    grad_L = [grad_Lx; grad_Ly];
    x_k = x_k - 0.01 * grad_L;  % Gradient descent update
    
    % Update Lagrange multipliers
    lambda = lambda + rho * h;
    mu = max(0, mu + rho * g);
    
    % Update penalty parameter if necessary
    if norm(h) > tol || any(g > tol)
        rho = 2 * rho;
    end
    
    % Check stopping criteria
    if norm(h) < tol && all(g <= tol)
        disp('Converged');
        break;
    end
end

% Output results
disp('Optimal x and y:');
disp(x_k);
