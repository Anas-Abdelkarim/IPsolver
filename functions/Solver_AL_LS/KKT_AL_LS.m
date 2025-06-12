function [KKT] = KKT_AL_LS(decision_variables,f_0,f_i,equality,parameters,option)
%KKT_LS(decision_variables,f,g,h,parameters,option);
% AL= Augumented Lagrangian; LS = least squares
%This function is run offline. All functions required for the online
%optimization is saved.
% decision_variables : the decision variables of the optimization  problem
% f_0 : is struct of three elements
% f_0.cost : the cost function
% f_0.error : cell of the error terms
% f_0.omega = the infromation matrix
% example J(X) = l(X) + e1' Omega1 e1 + e2' Omega2 e2
%  ....>  f_0.cost = l(X); f_0.error= {e1,e2}; f_0.omega = {Omega1,Omega2}
% if we have only l(X) we can say f_0 = l(X)
% f_i is matrix of  the inquality constraints
% equality is  matrix the inquality constraints
% decision_variables : the decision variables of the optimization  problem
% parameters         : the constants of the optimization problem
%%
% Check the options

slack_flag = 1; % true = Bazzana et al. ; flase = Sodhi et al. 
%% variables decleration
n =  length(decision_variables) ; % n is the number of the decision variablesq
q =  length(f_i);  % 	 is the number of the inqaulity constraints
l =  length(equality); % l is the number of the eqautilty constraints


lambda = sym('lambda',[q 1],'real');  % the dual variables for the inequaltiy constriants
active  = sym('active',[q 1],'real');   % 1 if g(x)>0  0 if g(x)<0 or 
ro_eq  = sym('ro_eq',[l 1],'real');     % the panelty for the equality constraints
ro_ineq = sym('ro_ineq',[q 1],'real');  % tthe panelty for the inequality constraints
gamma  = sym('gamma',[l 1],'real');     % the dual variables for the equaltiy constriants


%% the Newten step
if isrow(decision_variables)
    decision_variables = decision_variables';
end
X_ = decision_variables;
KKT_matrix = zeros(n,n);
KKT_vector = zeros(n,1);


%% * First and (second) derivates*

%1- cost function
if isstruct(f_0)
    if isfield(f_0, 'cost')
        f_0_cost = f_0.cost;
    else
        f_0_cost = 0;
    end
    error = f_0.error;
    information_matrix = f_0.omega;
else
    f_0_cost = f_0;
    information_matrix = [];

end


L = 0;
grad_L = 0;


if ~isempty(f_0_cost)
    jac_f_0_cost   = jacobian(f_0_cost,X_); % the first derivate of cost, i.e nebla_f_0
    hess_f_0_cost   = hessian(f_0_cost,X_)  ;
    KKT_matrix = KKT_matrix + .5*hess_f_0_cost ; % we multiply with .5 beacuse in FG they always devide by 2 the linear system
    KKT_vector = KKT_vector + -.5*jac_f_0_cost' ;
    f_0_tatal = 0 + f_0_cost;
end

%2 error terms (cost) (least squares)
for i = 1 : max(size(information_matrix))
    jac_error_cost = jacobian(error{1,i},X_);
    KKT_matrix = KKT_matrix + jac_error_cost'*information_matrix{1,i}*jac_error_cost;
    KKT_vector = KKT_vector + -jac_error_cost'*information_matrix{1,i}*error{1,i};
    f_0_tatal = f_0_tatal + error{1,i}'*information_matrix{1,i}*error{1,i};
end
L = f_0_tatal;
grad_L = jacobian(f_0_tatal,X_)';  % alway grad_L = -2 of KKT vector

% inequality constraints
if q>0
    f_i = lhs(f_i) - rhs(f_i);  % the inequality constraints in form f<=0
end
switch slack_flag
    case 1
        for i= 1 : q
            jac_error_ineq= jacobian(f_i(i),X_); % active if g(x)>-lambda/2*ro_ineqality
            omega_ineq = active(i)*ro_ineq(i);
            KKT_matrix = KKT_matrix + jac_error_ineq'*omega_ineq*jac_error_ineq;
            Weighted_error = -active(i)*(omega_ineq*f_i(i)+ .5*lambda(i));
            KKT_vector = KKT_vector + jac_error_ineq'*Weighted_error;
            f_i_p = active(i)*f_i(i) + (1-active(i))*(-.5*lambda(i)/ro_ineq(i));
            L = L +  omega_ineq * f_i_p^2 + lambda(i) * f_i_p  ;
            grad_L = grad_L + (-2)*jac_error_ineq'*Weighted_error;
        end

    case 0
        for i= 1 : q
            jac_error_ineq= jacobian(f_i(i),X_); 
            omega_ineq = active(i)*ro_ineq(i);   %active if g(x)>0
            KKT_matrix = KKT_matrix + jac_error_ineq'*omega_ineq*jac_error_ineq;
            Weighted_error = - (omega_ineq*f_i(i) + .5*lambda(i));
            KKT_vector = KKT_vector + jac_error_ineq'*Weighted_error;
            L = L +  omega_ineq * f_i(i)^2 + lambda(i) * f_i(i) ;
            grad_L = grad_L + (-2)*jac_error_ineq'*Weighted_error;
        end

end


% equality cosntraints

if l>0
    equality =  lhs(equality) - rhs(equality) ;  % the equality constraints in form g==0
end
for i= 1 : l
    jac_error_eq = jacobian(equality(i),X_);
    omega_eq = ro_eq(i);
    KKT_matrix = KKT_matrix + jac_error_eq'*omega_eq*jac_error_eq;
    Weighted_error = - (omega_eq*equality(i)+ 0.5*gamma(i));
    KKT_vector = KKT_vector + jac_error_eq'*Weighted_error;
    L = L + omega_eq * equality(i)^2 + gamma(i) * equality(i) ;
    grad_L = grad_L + (-2)*jac_error_eq'*Weighted_error;
end
 

%% Export as matlab functions
input= [decision_variables;parameters; active; lambda; gamma; ro_ineq; ro_eq];

KKT_matrix_func= sym2func(KKT_matrix,'Vars',input,option) ;
KKT_vector_func= sym2func(KKT_vector,'Vars',input,option);

f_0_func    = sym2func(f_0_tatal,'Vars',input,option); % for evaluation after finishing the calculations

grad_L_func  = matlabFunction(grad_L,'Vars',input); % for evaluation after finishing the calculations

if   ~isempty(f_i)
    f_i_func    = sym2func(f_i,'Vars',[decision_variables;parameters],option); % for evaluation after finishing the calculations
else
    f_i_func=[];
end

if   ~isempty(equality)
    equality_func    = sym2func(equality,'Vars',[decision_variables;parameters],option); % for evaluation after finishing the calculations
else
    equality_func=[];
end
 
%% return
KKT.KKT_matrix_func            = KKT_matrix_func;
KKT.KKT_vector_func            = KKT_vector_func;
KKT.f_0                        = f_0         ;% symbolic varibales for the  cost function
KKT.grad_L_func                = grad_L_func;
KKT.equality                   = equality ;% symbolic variables for the equality==0
KKT.f_0_func                   = f_0_func    ;
KKT.f_i_func                   = f_i_func    ;% function of the inquality (without slack variables)
KKT.input                      = input    ;
KKT.equality_func              = equality_func;
KKT.decision_variables         = decision_variables;
KKT.parameters                 = parameters        ;
KKT.option                     = option;
KKT.algorithm                  = 'PD_standard_LS';
KKT.facts                      = [n q l slack_flag]; %number of the decsion variable;inqualities and equalities
end