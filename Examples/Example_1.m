%%%%%%%%%%%%%%%%%%%       Example 1      %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              How to start with IPsolver 
%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We consider the optimization problem
% 
%     minimize     (x1-.5)^4 + 5*(x2-1)^4
%    subject to        x2     <=  14-2*x1
%                     .5*x1   <= -1.5+x2
%                      -1     <=  x1^2
%                      -1     <=  x2^2
%                     x2-3*x1  =  -6, 
% where a, b and c are parameters that are constants with respect to the  
% optimization problem, but we have the freedom to change them whenever 
% we call the solver.

clc
clear 

syms x1  x2     real
decision_variables= [x1 x2]  % the decision variables names of the optimization problem 

f_0     =    (x1-.5)^4 + 5*(x2-1)^4; % the cost function 
f_i     = [      x2       <=  14-2*x1
              .5*x1       <= -1.5+x2
                -1        <=  x1^2
                -1        <=  x2^2     ]; %the inequality constraints      
equality=[    x2-3*x1     ==  -6        ];  %the equality constraints.
%% call the solver
% Type the commad (help IPsolver) for more details

x_initial=[3 ,3];


algorithm = 'primal_dual_standard'   ;      % slack barrier solver: this algorithm does NOT require feasible point 
%algorithm    = 'primal_dual';          % primal dual with: this algorithm does NOT require feasible start point 

solution= IPsolver(x_initial,decision_variables,f_0,f_i,equality,algorithm);

% Extract the results
number_of_iterations       =  solution.num_iteration
x_optimal_solution         =  solution.x_optimal
minimum_cost               =  solution.cost_value
solver_time                =  solution.solver_time