function [KKT] = KKT(varargin);
% KKT(decision_variables,f_0,f_i,equality,parameter,algorithm,option)
% This function is called once for the optimization problem. It prepare the 
% Newton step.
% 1) decision variables are the decision variables in the optimization problem 
% defined by symbolic variables. Example: syms x y ; decision_variables=[x;y];
% 2)f_0 is the cost function. Example: f_0= 5*x^2 +2*y^2;
% 3)f_i is the inequality constraints. Example: f_i= [x+y<=100; 1<=x];
% 4)equality constraints of the optimization problem. Example
% equality=[x+y==0;x==3]; 
% 5)algorithm: leave it empty, we choose the defualt or choose one of the follows
% 1- slack barrier algorithm: 'slack_barrier'
% 2- primal dual with slack variables algorithm:'primal_dual'
% 2- primal dual algorithm :'primal_dual_standard' 
% 4- barrier algorithm: 'barrier' 
% Examples: algorithm= 'slack_barrier', algorithm= 'primal_dual', or
% algorithm= 'primal_dual_standard' or algorithm= 'barrier'.  
% 6)you can add some constants to optimization problem for example you
% can define f_0 = a*x^2 +b*y^2; and a or b is not decision variables 
% The parameters might be any number. They must be vector in vector form. Use reshape to 
% to convert the matrix to vector(vector_matrix =reshape(Your_matrix ,[],1))
% To define the parameters, do the follows: 
% syms a b; parameters= [ a; b]; and parameter_subs = [5;2]; 
% 8)option: type (help option_check) to know all details  
            
switch nargin
    case 0
         error('Not enough input arguments')
    case 1
         error('Not enough input arguments')
    case 2
        decision_variables =varargin{1}; 
        f_0                =varargin{2};
        f_i                =[];
        equality           =[];
        parameters         =[];       
        algorithm           =[];
        option             =[];
        disp('The default solver type that has been chosen is slack barrier (sb)') 
       
   case 3
        decision_variables  =varargin{1}; 
        f_0                 =varargin{2};
        f_i                 =varargin{3};
        equality            =[];
        parameters          =[];        
        algorithm           =[];
        option              =[];
        disp('The default solver type that has been chosen is slack barrier (sb)') 
   case 4
        decision_variables  =varargin{1}; 
        f_0                 =varargin{2};
        f_i                 =varargin{3};
        equality            =varargin{4};
        parameters          =[];        
        algorithm           =[];
        option              =[];
        disp('The default solver type that has been chosen is slack barrier (sb)') 
    case 5
        decision_variables  =varargin{1}; 
        f_0                 =varargin{2};
        f_i                 =varargin{3};
        equality            =varargin{4};           
        parameters          =varargin{5};   
        algorithm           =[];
        option              =[];
        disp('The default solver type that has been chosen is slack barrier (sb)') 
    case 6
        decision_variables  =varargin{1}; 
        f_0                 =varargin{2};
        f_i                 =varargin{3};
        equality            =varargin{4};
        parameters          =varargin{5}; 
        algorithm           =varargin{6};
        option              =[];
        disp('The default solver type that has been chosen is slack barrier (sb)') 
        case 7
        decision_variables  =varargin{1}; 
        f_0                 =varargin{2};
        f_i                 =varargin{3};
        equality            =varargin{4};
        parameters          =varargin{5};         
        algorithm           =varargin{6};
        option              =varargin{7} ;
    otherwise
        error('Too many input arguments')
end
%% Choose the default algorithm
if isempty(algorithm)
      algorithm ='primal_dual';
end

algorithm = algorithm_check(algorithm);


%% checking the options
 option=options_check(option,decision_variables);
  
 %%  converting the inputs to column vectors
assume(decision_variables, 'real')

if isrow(decision_variables)
    decision_variables =sym2cell(decision_variables)';
    decision_variables=cell2sym(decision_variables);
end

if ~isempty(parameters)
   assume(parameters, 'real')
if isrow(parameters) 
   parameters =sym2cell(parameters)';
    parameters=cell2sym(parameters);
  end
end


if isrow(f_i)
    f_i =sym2cell(f_i)';
    f_i=cell2sym(f_i);

end



if isrow(equality) 
   equality =sym2cell(equality)';
    equality=cell2sym(equality);
end


%%
switch algorithm
    case 'slack_barrier'
        [KKT] = KKT_slack_barrier(decision_variables,f_0,f_i,equality,parameters,option);
    case 'primal_dual'
        [KKT] = KKT_primal_dual(decision_variables,f_0,f_i,equality,parameters,option);
    case 'primal_dual_standard'
        [KKT] = KKT_primal_dual_standard(decision_variables,f_0,f_i,equality,parameters,option);  
    case 'PD_standard_LS'
        [KKT] = KKT_PD_standard_LS(decision_variables,f_0,f_i,equality,parameters,option);  
    case 'barrier'
        [KKT] = KKT_barrier(decision_variables,f_0,f_i,equality,parameters,option);
      
    otherwise
      error("algorithm must be one of the following: 'PD_standard_LS', 'primal_dual','slack_barrier', 'barrier', or 'slack_barrier'. Type (help IPsolver) for more details")
  
end

KKT.algorithm           = algorithm; 

end

 




