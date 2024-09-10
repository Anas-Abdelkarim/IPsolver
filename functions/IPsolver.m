    function [IPsolver] = IPsolver(varargin);
   % To call this function there is two methods depends on the situation
    %%%%%%%%%%%%%%%%%%%% option 1 %%%%%%%%%%%%%%%%%%%%%%%%%%
    % IPsolver(x_initial,decision_variables,f_0,f_i,eqaulity,algorithm ,option,parameters,parameters_subs,warm_start);
    % this options is used if you want the solver to call KKT function
    % automatically. KKT function is explained in KKT.m
    %%%%%%%%%%%%%%%%%%% option 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IPsolver(KKT1,x_initial,parameters_subs,warm_start)
    % This option is used when the KKT function is called offline (we need
    % to call it only once to prepair the the Newton step). Feed the KKT to
    % this function as first argument. 
    % This option is helpful in Model  predictive control problems. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       INPUT ARGUMENTS                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) decision variables are the decision variables in the optimiazion problem 
    % based on symbolic variables. Example: syms x y ; decision_variables=[x;y];
    % 2)x_initial is vector contains the initial values of the variables
    % ex: x_initial = [ 1 ;2]; 
    % 3)f_0 is the cost function. Example: f_0= 5*x^2 +2*y^2;
    % 4)f_i is the inequality constraints. Example: f_i= [x+y<=100; 1<=x];
    % 5)equality constraints of the optimization problem. Example
    % equality=[x+y==0;x==3]; 
    % 6)algorithm: leave it empty, we choose the defualt or choose one of the follows
    % 1- slack barrier algorithm: 'slack_barrier'
    % 2- primal dual with slack variables algorithm:'primal_dual'
    % 2- primal dual algorithm :'primal_dual_standard' 
    % 4- barrier algorithm: 'barrier' 
    % Examples: algorithm= 'slack_barrier', algorithm= 'primal_dual', or
    % algorithm= 'primal_dual_standard' or algorithm= 'barrier'.  
    % 7)you can add some constants to optimization problem for example you
    % can define f_0 = a*x^2 +b*y^2; and a or b is not decision variables 
    % to convert the matrix to vector(vector_matrix =reshape(Your_matrix ,[],1))
    % To define the parameters, do the follows: 
    % syms a b; parameters= [ a; b]; and parameter_subs = [5;2]; 
    % 8)option: type (help option_check) to know all details  
    %8)warm_start to have the optimal values of the last solution of the optimization 
    % problem. Following choices are valid for warm_start depending on the
    % solver type: warm_start.x_optimal, warm_start.lambda_optimal,
    % warm_start.mu_optimal, warm_start.gamma; warm_start.slack    
    % 10)KKT1: type (help KKT) to know all details. In general call it using: 
    % KKT(decision_variables,f_0,f_i,equality,parameters,algorithm       ,option)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       OUTPUT ARGUMENTS                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1)IPsolver.x_optimal     : the solution of the optimization problem
    %2)IPsolver.gamma_optimal : the optimal equality dual variables
    %3)IPsolver.cost_value    : the function value
    %4)IPsolver.num_iteration : number of the iterations to find the
    %optimala value
    %5)IPsolver.s_record      : the step size at each iteration
    %6)IPsolver.Newton_step_record : Newton step at each iteration
    %7)IPsolver.t_record      : t is the approximation of the barrier
    %function. At t increase we approach the ideal barrier function
    %8)IPsolver.x_record      :  the value of x at each iteration
    %9)IPsolver.gamma_record  : the value of equality dual variables at each iteration
    %10)IPsolver.lambda_record: the value of inequality dual variables at each iteration
    %11)IPsolver.slack_record : the value of slack variables at each iteration
    %12)IPsolver.slack_optimal: the optimal value of slack  variables 
    %13)IPsolver.lambda_optimal: the optimal inequality dual variables
    %14)IPsolver.mu_optimal    : the optimal equality dual variables (for
    %equality constraints with slack variables)
    %15)IPsolver.search_x_start: detail about finding feasible x initial
    %16)IPsolver.solver_time : the calculation time  required (not
    %including building the KKT system because this can be done once and
    %separately)

    
    
    
if isa(varargin{1},'double')
    switch nargin
        case 1
               error('Not enough input arguments')
        case 2
               error('Not enough input arguments')
        case 3 
            x_initial          =varargin{1};
            decision_variables =varargin{2};
            f_0                =varargin{3};
            f_i                =[];
            equality           =[];
            algorithm          =[];
            option             =[];
            parameters         =[];
            parameters_subs    =[];
            warm_start         =[];
        case 4
            x_initial          =varargin{1};
            decision_variables =varargin{2};
            f_0                =varargin{3};
            f_i                =varargin{4};
            equality           =[];
            algorithm          =[];
            option             =[];
            parameters         =[];
            parameters_subs    =[];
            warm_start         =[];
        case 5
            x_initial          =varargin{1};
            decision_variables =varargin{2};
            f_0                =varargin{3};
            f_i                =varargin{4};
            equality           =varargin{5};
            algorithm          =[];
            option             =[];
            parameters         =[];
            parameters_subs    =[];
            warm_start=[];
        case 6
            x_initial            =varargin{1};
            decision_variables   =varargin{2};
            f_0                  =varargin{3};
            f_i                  =varargin{4};
            equality             =varargin{5};
            algorithm            =varargin{6};
            option               =[];
            parameters           =[];
            parameters_subs      =[];
            warm_start           =[];
        case 7
            x_initial            =varargin{1};
            decision_variables   =varargin{2};
            f_0                  =varargin{3};
            f_i                  =varargin{4};
            equality             =varargin{5};
            algorithm            =varargin{6};
            option               =varargin{7};
            parameters           =[];
            parameters_subs      =[];
            warm_start           =[];
        case 8
                error('Not enough input arguments. check parameters_subs')
        case 9
            x_initial          =varargin{1};
            decision_variables =varargin{2};
            f_0                =varargin{3};
            f_i                =varargin{4};
            equality           =varargin{5};
            algorithm          =varargin{6};
            option             =varargin{7};
            parameters         =varargin{8};
            parameters_subs    =varargin{9};
            warm_start         =[];
        case 10
            x_initial          =varargin{1};
            decision_variables =varargin{2};
            f_0                =varargin{3};
            f_i                =varargin{4};
            equality           =varargin{5};
            algorithm          =varargin{6};
            option             =varargin{7};
            parameters         =varargin{8};
            parameters_subs    =varargin{9};
            warm_start         =varargin{10};   
    end
  
   if  isempty(algorithm)       
       algorithm='primal_dual';
   end      
 algorithm = algorithm_check(algorithm);
  
  if ~isempty(parameters)
      f_0= subs(f_0,parameters,parameters_subs)  ;
      f_i= subs(f_i,parameters,parameters_subs)  ;
      equality= subs(equality,parameters,parameters_subs)  ;
      parameters=[];
      parameters_subs=[];
  end
  
  [KKT1] = KKT(decision_variables,f_0,f_i,equality,parameters,algorithm,option);   


  
    
else
    
switch nargin
    case 0
        error('Not enough input arguments')
    case 1
        KKT1             =varargin{1};
        x_initial       =[];
        parameters_subs =[];
        warm_start=[];
    case 2
        KKT1             =varargin{1};
        x_initial       =varargin{2};
        parameters_subs =[];
        warm_start=[];
    case 3
        KKT1             =varargin{1};
        x_initial       =varargin{2};
        parameters_subs =varargin{3};
        warm_start=[];
    case 4
      KKT1             =varargin{1};
      x_initial       =varargin{2};
      parameters_subs =varargin{3};
      warm_start=varargin{4};
    otherwise
        error('Too many input arguments')
end
end

if isrow(x_initial) 
    x_initial=x_initial';
end

if ~isempty(x_initial)
if length(x_initial)~=length(KKT1.decision_variables)
    error('Number the decision variables must be equale to length x_initial');
end
end
if isrow(parameters_subs) 
    parameters_subs=parameters_subs';
end

switch KKT1.algorithm     
    case 'slack_barrier'
        [IPsolver] = solver_slack_barrier(KKT1,x_initial,parameters_subs,warm_start);
    case 'primal_dual'    
        [IPsolver] = solver_primal_dual(KKT1,x_initial,parameters_subs,warm_start);    
     case 'primal_dual_standard'    
        [IPsolver] = solver_primal_dual_standard(KKT1,x_initial,parameters_subs,warm_start);
     case 'PD_standard_LS'    
        [IPsolver] = solver_PD_standard_LS(KKT1,x_initial,parameters_subs,warm_start);
     case 'AL_LS'    
        [IPsolver] = solver_AL_LS(KKT1,x_initial,parameters_subs,warm_start);
     case 'barrier'
        [IPsolver] = solver_barrier(KKT1,x_initial,parameters_subs,warm_start);
    otherwise
      error("algorithm must be one of the following:'AL_LS', 'PD_standard_LS', 'primal_dual','slack_barrier', 'barrier', or 'slack_barrier'. Type (help IPsolver) for more details")

end
 end
    