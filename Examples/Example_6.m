%%%%%%%%%%%%%%%%%%%       Example 1      %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              How to start with IPsolver
%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We consider the optimization problem
%
%     minimize     (x-2)^2 +  (y-9)^2 + (z-50)^2
%    subject to
%                      y+z    == 3,
%                      x+y    == 2.5,
%                      x      <= 1
%                      y      <= 4
%                      z      <= 5
% 
% where a, b and c are parameters that are constants with respect to the
% optimization problem, but we have the freedom to change them whenever
% we call the solver.

clc
clear

syms x  y  z nu nu2 s1 s2  k1 k2 real
flag_LS_cost = 0;
decision_variables= [x y z];  % the decision variables names of the optimization problem

x_initial = [-40, -20, -10];


    switch flag_LS_cost
        case 0
            f_0     =    (x-2)^2 + (y-9)^2 + (z-50)^2; % the cost function
        case 1
            decision_variables = [decision_variables nu nu2]
            x_initial = [x_initial, 0, 0 ];
            f_0.error = {x-2,y-9,z-50, [y+z-3 ;x+y-2.5 ; nu ;  nu2] };
            f_0.omega = {1,1,1, [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0]};
    end

 
f_i     =  [x <=  1 ; y <= 4;    z  <=5];
%f_i = [x<=05];
% f_i = [x <=  5]
  
equality= [y+z    == 3;   x+y    == 2.5];           
%equality = [y+z    == 3];
  

%% call the solver
% Type the commad (help IPsolver) for more details

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
solution.lambda_record
solution.x_record
solution.gamma_optimal
solution.x_optimal
size(solution.x_record)

%solution.gamma_record



%{
solution.lambda_optimal
solution.gamma_optimal
solution.x_optimal

%}
%num_iteration_Newton       = solution.num_iteration_Newton