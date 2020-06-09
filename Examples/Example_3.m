%%%%%%%%%%%%%%%%%%%       Example 3      %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example shows how to deal with parameters in matrix form. 
% Moreover, this example presents the available algorithms in the solver.
%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We consider the optimization problem
% 
%     minimize     (x(1)-.5)^4 + u*(x(2)-1)^4
%    subject to         -1   <=x(1)^2
%                       -1   <=x(2)^2
%                      A*x   <=  b
%                   
%                     
% where  x=[x(1)        A= [-3  1       b = [-6         u=5
%           x(2)];           2  1             14
%                           .5 -1];          -1.5]


clc  
clear  

x = sym('x', [2 1], 'real');

A = sym('A', [3 2], 'real');
b = sym('b', [3 1], 'real');
syms u real; 

decision_variables= x;  % the decision variables names of the optimization problem 
parameters        = [reshape(A,[],1);b;u];  % to convert the matrix to vector form we use reshape

f_0 =  (x(1)-.5)^4 + u*(x(2)-1)^4; % the cost function 
f_i=[    -1        <=x(1)^2
         -1        <=x(2)^2
         A*x       <=  b     ]; %the inequality constraints  %the inequality constraints

equality=[];       %the equality constraints is empty

% algorithm = 'slack_barrier'   ;      % slack barrier solver: this algorithm does NOT require feasible point 
%algorithm = 'primal_dual_standard' ; % primal dual solver: this algorithm does require feasible start  point 
algorithm = 'barrier'  ;             % barrier solver: this algorithm does require feasible start point 
% algorithm    = 'primal_dual';          % primal dual with: this algorithm does NOT require feasible start point 
%% 1- Newton step :
  Option.elimination=1; % this option builds the Newton step based on the reduced KKT system (save time)
  Option.data_recording=1; % with this option the solver does not record ever iteration data (save time and memory)
     
  NewtonStep = KKT(decision_variables,f_0,f_i,equality,parameters,algorithm,Option);

  %% call the solver
  % online calculation

x_initial =[4;5];    
A_subs= [-3  1 
          2  1 
         .5 -1];
b_subs =[-6 ;+14;-1.5]; 
u_subs =   5  ;

parameters_subs = [reshape(A_subs,[],1) ;b_subs;u_subs]; 
solver_IPsolver = IPsolver(NewtonStep,x_initial,parameters_subs); 

%% 1- Comparing with fmincon:

% we have built an interface for fmincon with symbolic variables    
 [solver_fmincon] =  fmincon_interface(decision_variables,f_0, f_i,equality,x_initial,parameters,parameters_subs);    

 %% Results harvesting 
% 1-  IPsovler 
n_iterations   =  solver_IPsolver.num_iteration;
x_optimal     =  solver_IPsolver.x_optimal;
cost_value     =  solver_IPsolver.cost_value;
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('1- IPsolver')
disp(["Number of the iterations:"+n_iterations])
text_show= "The solutoin of the optimization problem: x1="+ x_optimal(1) +" and x2=" +x_optimal(2) ;
disp(text_show);
disp(["The cost value:"+cost_value])
text_show = "The calculation time for is: ("+solver_IPsolver.solver_time*1000 +") ms";
disp(text_show);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
% 2- fmincon solver using interior-point algorithm
n_iterations_fmincon  =  solver_fmincon.num_iteration;
x_optimal_fmincon     =  solver_fmincon.x_optimal;
cost_value_fmincon    =  solver_fmincon.cost_value;
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('2- fmincon solver using interior-point algorithm')
disp(["Number of the iterations:"+n_iterations_fmincon])
text_show= "The solutoin of the optimization problem: x1="+ x_optimal_fmincon(1) +" and x2=" +x_optimal_fmincon(2) ;
disp(text_show);
disp(["The cost value:"+cost_value_fmincon])
text_show = "The calculation time is ("+solver_fmincon.solver_time*1000 +") ms";
disp(text_show);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

 
% constraints drawing
f=f_0;
f=subs(f_0,parameters,parameters_subs);

p1=fcontour((f),[-1,6.5],'LineWidth', 1);
levels=logspace((log10(0.1)),(log10(5500)),15)';
p1.LevelList = levels;
hold on
scat1=scatter([3,4,5],[3,6,4],5,'g','filled');
scat3=scatter([.5],[1],'bo','filled');

t1= text(3.09,3.01,'(3,3)','FontSize',10);
t1= text(4.05,6.04,'(4,6)','FontSize',10);
t1= text(5.05,4.04,'(5,4)','FontSize',10);


%t1= text(0.55,1.04,'(0.5,1)','FontSize',10);

A = [3,4,5] ;% x coordinates ;
B=  [3,6,4] ;% y coordintes;
C = [ -13,13] ;% x coordinates ;
%l3= plot(A,B,'black' ,'LineWidth', 1);
polytope= patch(A,B,'g');
polytope.EdgeColor='none';
h = legend(p1 ,'Location', 'Best');
 h.Interpreter = 'tex';
 axis equal
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
 % set(gcf, 'Units', 'Normalized', 'Position', [0, 0.04, .7, .6]);
bar= colorbar;
colormap(jet)
plot=gca;
plot.XLim = [-1 6.5];
plot.YLim = [-1 6.5];
 plot.Units='Pixels';
 shg
 
%}