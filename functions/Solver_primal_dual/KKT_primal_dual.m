function [KKT] = KKT_primal_dual(decision_variables,f,g,h,parameters,option);
%KKT_sb(decision_variables,f,g,h,parameters,option);.
%This function is run offline. All functions required for the online
%optimization is saved.
% decision_variables : the decision variables of the optimization problem
% f is the cost function
% g is  the inquality constraints 
% h is  the inquality constraints 
% decision_variables : the decision variables of the optimization problem
% parameters         : the constants of the optimization problem.
%% variables decleration


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

n =  length(decision_variables) ; % n is the number of the decision variablesq
q =  length(g);  % 	 is the number of the inqaulity constraints 
l =  length(h); % l is the number of the eqautilty constraints
slack = sym('slack',[q 1],'real');  % the options_check variables for the inequaltiy constriants
lambda    = sym('lambda',[q 1],'real');  % the dual variables for the equaltiy constriants with slack variables
gamma = sym('gamma',[l 1],'real');     % the dual variables for the equaltiy 
syms t       real; 
Delta_x      = sym('Delta_x',[n 1],'real');  % the update values of the decision variables
Delta_slack  = sym('Delta_slack',[q 1],'real');  % the update values of the decision variables
Delta_lambda     = sym('Delta_lambda',[q 1],'real');  % the update values of the decision variables
Delta_gamma   = sym('Delta_gamma',[l 1],'real') ; % the update values of the decision variables
%% * First and (second) derivates*
%cost function
grad_f   = jacobian(f,decision_variables)' ; % the first derivate of cost, i.e nebla_f
hess_f   =  hessian(f,decision_variables)  ; % the second derivate of the cost funciton

%inequality constraints
if   ~isempty(g) 
g = lhs(g) - rhs(g);  % the inequality constraints become equality with salck variables
Dg       = jacobian(g,decision_variables) ; % the first derivate of inequality with respect to x 

   if ~isempty(symvar(Dg))  %means some inequality has symbolic variables
      assume(Dg,'real') 
   end 
   
%%%% to test if all inequalities are linear
IneqLinAll= isempty(symvar(subs(Dg,parameters,ones(length(parameters),1)))); % 1 means all linear inequality
%%%%%

end
%equality constraints
if   ~isempty(h)
    h    =   lhs(h) - rhs(h); %the equaility constraints in form Ax-b=0 
    Dh = jacobian(h,decision_variables) ; % the first derivate of equaility, i.e A matrix
    if ~isempty(symvar(Dh)) % means some equality has symbolic variables
        assume(Dh,'real')
    end
         
else
    Dh=[];
end
%% the Newten step

daulHess_g = zeros(n,n);
if   ~isempty(g) & ~IneqLinAll
  for i = 1:q
    daulHess_g =   daulHess_g +lambda(i)*hessian(g(i),decision_variables);
   end
end
H_1 = hess_f+daulHess_g;
if   ~isempty(g)
     if elimination_flag == 0
        
     KKT_matrix =[    H_1      ,    Dg'     ,   Dh'       , zeros(n,q)
                    zeros(q,n) , diag(slack),  zeros(q,l) ,diag(lambda)
                      Dh       ,  zeros(l,q), zeros(l,l)  , zeros(l,q)
                      Dg       ,  zeros(q,q) , zeros(q,l)   , eye(q,q)   ];
                  
                  
     Delta_slack  = [] ;
     Delta_lambda_step= [];
     else  %% elimination_flag == 1
          A1= Dg'*diag(slack.^-1)*diag(lambda);
          A2= A1*(g+1/t*(lambda.^-1));
          KKT_matrix =[    H_1+A1*Dg         ,    Dh'    
                               Dh            , zeros(l,l) ];       
         
  
     Delta_lambda_step= diag(slack.^-1)*diag(lambda)*(g+1/t*(lambda.^-1)+Dg*Delta_x);
     Delta_slack = -diag(lambda.^-1)*(diag(slack)*Delta_lambda+diag(lambda)*slack-1/t ) ;     
     end 
                    
    
    
else
     KKT_matrix =[ hess_f   ,    Dh'
                    Dh    , zeros(l,l) ];                             
     Delta_slack  = [] ;
     Delta_lambda_step= [];
end
% KKT vector
%1)###### r_dual
r_dual      =  grad_f ;
if   ~isempty(g) 
  r_dual     = r_dual+Dg'*lambda ;
end
if   ~isempty(h) 
  r_dual   =  r_dual +Dh'*gamma ; 
end


%2)###### r_cent
if   ~isempty(g) 
   r_cent   = diag(lambda)*slack-(1/t)       ;
else
    r_cent = [];
end

%3)###### r_h
if   ~isempty(h) 
   r_pri = h    ;
else
    r_pri =[];
end

%3)###### r_g
if   ~isempty(g)
    r_g = g+slack;
   else 
   r_g =[];
end


r_t    = [r_dual;r_cent;r_pri;r_g];


if elimination_flag == 0 | isempty(g)
        KKT_vector= - r_t;     
end

if elimination_flag == 1
    if ~isempty(g)       
        KKT_vector= -[r_dual+A2;r_pri];          
    end
end
         
 


%% Preparation for prime daul interior point solver
input= [decision_variables;parameters;slack;lambda;gamma;t];
KKT_matrix_func= sym2func(KKT_matrix,'Vars',input,option) ;
KKT_vector_func= sym2func(KKT_vector,'Vars',input,option);

    
% residuals vectors
r_t_func     = sym2func(r_t,'Vars',input,option);
r_dual_func  = sym2func(r_dual,'Vars',input,option);
r_pri_func   = sym2func([r_pri;r_g],'Vars',input,option);
% The handle function for cost funtion equality and dynamic cost function 
f_func    = sym2func(f,'Vars',input,option); % for evaluation after finishing the calculations

if elimination_flag ==1;   
        inputForDelta_lambda= [input;Delta_x];
        Delta_lambda_func= sym2func(Delta_lambda_step,'Vars',inputForDelta_lambda,option);   
        inputForDelta_slack= [input;Delta_lambda];      
        Delta_slack_func= sym2func(Delta_slack,'Vars',inputForDelta_slack,option);
else
       Delta_lambda_func= sym2func([],'Vars',input,option);  
       Delta_slack_func= sym2func([],'Vars',input,option);
end

if   ~isempty(g)
eta_hat = (slack)'*(lambda);
else 
eta_hat =1e-100 ;
end
eta_hat_func= sym2func(eta_hat,'Vars',input,option);
 g_func    = sym2func(g,'Vars',input,option); % for evaluation after finishing the calculations
% if   ~isempty(g) % means we have inequality constraints 
%   g_func    = sym2func(g,'Vars',input,option); % for evaluation after finishing the calculations
% else
%  g_func   =  [];
% end
%% return
KKT.KKT_matrix_func            = KKT_matrix_func;
KKT.KKT_vector_func            = KKT_vector_func;
KKT.Delta_slack_func           = Delta_slack_func;
KKT.Delta_lambda_func          = Delta_lambda_func   ;
KKT.eta_hat_func               = eta_hat_func;
KKT.r_t_func                   = r_t_func    ; 
KKT.r_dual_func                = r_dual_func ;
KKT.r_pri_func                 = r_pri_func  ; 
KKT.f                          = f           ;% cost function (symbolic varibales form)
KKT.f_i                        = g           ;% inequality function without conversion (symbolic varibales form)
KKT.equality                   = h           ;% equality function without conversion (symbolic varibales form)
KKT.f_func                     = f_func      ;
KKT.g_func                     = g_func      ; %inequality function of slack variables (symbolic varibales form)
KKT.option                     = option      ; %inequality function of slack variables (symbolic varibales form)
KKT.decision_variables         = decision_variables;
KKT.parameters                 = parameters;
KKT.elimination_flag           = elimination_flag   ;
KKT.algorithm                  ='primal_dual'; 
KKT.facts                      = [n q l]; %number of the decsion variable;inqualities and equalities
end