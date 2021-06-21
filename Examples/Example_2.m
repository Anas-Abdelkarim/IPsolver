%%%%%%%%%%%%%%%%%%%       Example 2      %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example shows how to prepare the Newton step and before calling the 
% solver. Also, it shows alternatives to initialize the decision variables
%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We consider the optimization problem
% 
%     minimize     (x1-.5)^4 + a*(x2-1)^4
%    subject to        x2      <=  14-2*x1
%                     c*x1     <= -1.5+x2
%                      -1      <=  x1^2
%                      -1      <=  x2^2
%                    x2-3*x1    =  b

 clc 
 clear  


syms x1  x2   a b c   real
decision_variables = [x1 x2]; % the decision variables names of the optimization problem 
parameters         = [a b c]; 

f_0      =    (x1-.5)^4 + a*(x2-1)^4     ; % the cost function 
f_i      = [     x2        <=  14-2*x1
                c*x1       <= -1.5+x2
                -1         <=  x1^2
                -1         <=  x2^2      ]; %the inequality constraints      
equality = [   x2-3*x1     ==  b         ];  %the equality constraints.
           
%% 1- Preparing the Newton step 
% Type the commad (help KKT) for more details 

Newten_step= KKT(decision_variables,f_0,f_i,equality,parameters);
%% 2- call the solver

% case 1: let the parameters same as in the Example_1 
% a = 5, b= -6 and c= .5

x_initial = [];
parameters_subs = [5 -6 .5];
case_1=  IPsolver(Newten_step,x_initial,parameters_subs);

% Extract the results
number_of_iterations_case1 =  case_1.num_iteration
x_optimal_case1            =  case_1.x_optimal
minimum_cost_case1         =  case_1.cost_value
solver_time_case1          =  case_1.solver_time


% case 2: let the parameters a = 8, b= -8 and c= .8

x_initial = x_optimal_case1; % here we use the optimal value obtained form the previous case
parameters_subs = [8 -3 .8];
case_2=  IPsolver(Newten_step,x_initial,parameters_subs);

% Extract the results
number_of_iterations_case2 =  case_2.num_iteration
x_optimal_case2            =  case_2.x_optimal
minimum_cost_case2         =  case_2.cost_value
solver_time_case2          =  case_2.solver_time


% case 3: let the parameters same as in the Example_1 
% a = 10, b= 1 and c= 0

x_initial = [];
warm_startr= case_2;  % Note in warm start, we provide initial values for not only the decision variables but also for the slack and dual variabls
parameters_subs = [2 -7 .6];
case_3=  IPsolver(Newten_step,x_initial,parameters_subs);

% Extract the results
number_of_iterations_case3 =  case_3.num_iteration
x_optimal_case3            =  case_3.x_optimal
minimum_cost_case3         =  case_3.cost_value
solver_time_case3          =  case_3.solver_time