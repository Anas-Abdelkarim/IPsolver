%%%%%%%%%%%%%%%%%%%       Example 4      %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%h%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%i%%%%%%%%%%%%%%%%%%%%%%

%We consider the optimization problem
% 
%     minimize          sqrt(x)
%    subject to         0    <=x
%                       x    <=10
%                
clc  
clear  

syms x;
decision_variables = x;  % the decision variables names of the optimization problem 
x_initial          = 1500;  % infeasible start point
 
f_0= -sqrt(x);
f_i= [0<=x;
      x<=10]; 
solution=IPsolver(x_initial,decision_variables,f_0,f_i);
% results %
solution.x_optimal
solution.cost_value


