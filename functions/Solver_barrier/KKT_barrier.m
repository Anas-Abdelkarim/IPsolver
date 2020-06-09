function [KKT] = KKT_barrier(decision_variables,f_0,f_i,equality,parameters,option);
%KKT_sb(decision_variables,f,g,h,parameters,option);.
%This function is run offline. All functions required for the online
%optimization is saved.
% decision_variables : the decision variables of the optimization problem
% f_0 is the cost function
% f_i is  the inquality constraints 
% equality is  the inquality constraints 
% decision_variables : the decision variables of the optimization problem
% parameters         : the constants of the optimization problem



for check_options = 1
if ~exist('option')
  option= options_check([]);
else
  if~isfield(option,'check_done')
     option= options_check(option);
   end
end
if option.check_done==0
    option= options_check(option);
end

if option.elimination 
   elimination_flag=1;
else 
   elimination_flag=0;
end


end

%% variables decleration
n =  length(decision_variables) ; % n is the number of the decision variables
q =  length(f_i);  % q is the number of the inqaulity constraints 
l =  length(equality); % p is the number of the eqautilty constraints
syms t       real; 
%%  First and (second) derivates*
%inequaility constraints
if   ~isempty(f_i)
f = lhs(f_i) - rhs(f_i) ;  % the inequaility constraints in form f<=0 
end
if   ~isempty(f)
f_0_bar = f_0-sum(log(-f)/t) ;% approximation of the inequalities 
end
grad_f_0_bar   = jacobian(f_0_bar,decision_variables)' ; % the first derivate of cost, i.e nebla_f
hess_f_0_bar   =  hessian(f_0_bar,decision_variables)  ;
%equaility constraints

if   ~isempty(equality)
    equality             = lhs(equality) - rhs(equality); %the equaility constraints in form Ax-b=0 
    A                    = jacobian(equality,decision_variables) ; % the first derivate of equaility, i.e A matrix.
    if ~isempty(symvar(A))  %means some equality has symbolic variables
       assume(A,'real') 
    end 
else
    A=[];
    equality=[];
end
%% the Newten step

     KKT_matrix = [ hess_f_0_bar         A'   
                            A     zeros(l,l)] ;
          
     KKT_vector = -[grad_f_0_bar; zeros(l,1)] ;
     
     Delta_x_bar      = sym('Delta_x_bar',[n 1],'real');  % the update values of the decision variables
     decrement=  Delta_x_bar'*grad_f_0_bar ;



%%%%% Finding feasible start point %%%%%%% 
KKT_x_start=[];
%% this done by preparing the KKT system of the OP 
if option.find_feasible_point 
      option.find_feasible_point=0; % to avoid infinite loop
    if ~isempty (equality)
        b=subs(equality,decision_variables,zeros(n,1));
        x_start_search.A= A;
        x_start_search.b= b;
      else
        x_start_search= [];
    end
   
    if ~isempty (f_i)
      syms level_unique real % 
      f_0_x_start= level_unique;
      f_i_x_start = f<= level_unique;
      equality_x_start= equality==0;
      decision_variables_x_start= [ level_unique;decision_variables];     
       KKT_x_start = KKT_barrier(decision_variables_x_start,f_0_x_start,f_i_x_start,equality_x_start,parameters,option);
     end
    
     if ~isempty(x_start_search)
      KKT_x_start.x_start_search.A_func= sym2func( x_start_search.A,'Vars',[decision_variables;parameters],option) ;
      KKT_x_start.x_start_search.b_func= sym2func( x_start_search.b,'Vars',[decision_variables;parameters],option) ;
     else
       KKT_x_start.x_start_search=x_start_search;
       KKT.KKT_x_start =KKT_x_start;
     end
      
       option.find_feasible_point=1; % to return it to its correct value

end
    
    
  
 
%% Preparation for the interior point solver

input= [decision_variables;parameters;t];
% the newton step
KKT_matrix_func= sym2func(KKT_matrix,'Vars',input,option) ;
%KKT_matrix_factorized.L= matlabFunction(L,'Vars',input,option);
%KKT_matrix_factorized.U= matlabFunction(U,'Vars',input,option);
%KKT_matrix_factorized.P= P                             ;
KKT_vector_func= sym2func(KKT_vector,'Vars',input,option);
% The cost function 
f_0_func    = sym2func(f_0,'Vars',input,option); % for evaluation after finishing the calculations
f_0_bar_func      = sym2func(f_0_bar,'Vars',input,option); % the approximated cost function  
f_i_func  = sym2func(f,'Vars',input,option); % for backtracking
equality_func    = sym2func(equality,'Vars',input,option); % 
inputForDecrement= [input;Delta_x_bar];
decrement_func= sym2func(decrement,'Vars',inputForDecrement,option);

  
%% return
KKT.KKT_matrix_func             = KKT_matrix_func;
%KKT.KKT_matrix_factorized      = KKT_matrix_factorized;
KKT.KKT_vector_func             = KKT_vector_func;
KKT.f_0_func                    = f_0_func   ; % for evaluation after finishing the calculations
KKT.f_i_func                    = f_i_func   ;  % for evaluation after finishing the calculations
KKT.f_0_bar_func                = f_0_bar_func   ;
KKT.equality_func               = equality_func  ;
KKT.decrement_func              = decrement_func ;
KKT.f_0                         = f_0            ;
KKT.f_i                         = f_i            ;
KKT.equality                    = equality       ;
KKT.decision_variables          = decision_variables;
KKT.parameters                  = parameters        ;
KKT.elimination_flag            = elimination_flag  ;
KKT.option                      = option         ; 
KKT.KKT_x_start                 = KKT_x_start    ;
KKT.algorithm                   = 'barrier'      ; 
KKT.facts                       = [n q l ]; %number of the decsion variable;inqualities and equalities
end