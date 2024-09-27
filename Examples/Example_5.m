%%%%%%%%%%%%%%%%%%%       Example 1      %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              How to start with IPsolver
%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We consider the optimization problem
%
%     minimize     (x-25)^2 + (y-50)^2
%    subject to

%                      x+y    =  20,
% where a, b and c are parameters that are constants with respect to the
% optimization problem, but we have the freedom to change them whenever
% we call the solver.

clc
clear

syms x  y   s1 s2  k1 k2 real
flag_LS_cost = 0;

switch flag_LS_cost
    case 0
        f_0     =    (x-50)^2 + (y-25)^2; % the cost function
    case 1
        f_0.error = {x-50,y-25};
        f_0.omega = {1,1};

end

 
f_i     = [ x <=  15
            y>=7
                  ];

equality=     x+y ==  20;                        %the equality constraints.
 

%% call the solver
% Type the commad (help IPsolver) for more details
decision_variables= [x y];  % the decision variables names of the optimization problem
x_initial=[-100,1200];

switch 1

    case 1
        algorithm = 'barrier_LS'; 
    case 2
        algorithm = 'PD_standard_LS';
    case 3
        algorithm = 'AL_LS';
    case 4
        algorithm = 'barrier';
    case 5
        algorithm = 'primal_dual_standard'   ;
    case 6
        algorithm = 'primal_dual'   ;
  
end

solution= IPsolver(x_initial,decision_variables,f_0,f_i,equality,algorithm);

% Extract the results
number_of_iterations       =  solution.num_iteration
x_optimal_solution         =  solution.x_optimal
minimum_cost               =  solution.cost_value
solver_time                =  solution.solver_time
%{
solution.lambda_optimal
solution.gamma_optimal
solution.x_optimal

%}
%num_iteration_Newton       = solution.num_iteration_Newton