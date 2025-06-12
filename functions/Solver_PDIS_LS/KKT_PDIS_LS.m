function [KKT] = KKT_PDIS_LS(decision_variables,f_0,f_i,equality,parameters,option)
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
        warning('Elimination flag in PDIS is always set to 1')
    end
end

n =  length(decision_variables) ; % n is the number of the decision variablesq
q =  length(f_i);  % 	 is the number of the inqaulity constraints
l =  length(equality); % l is the number of the eqautilty constraints
slack = sym('slack',[q 1],'real');  % the options_check variables for the inequaltiy constriants
lambda    = sym('lambda',[q 1],'real');  % the dual variables for the equaltiy constriants with slack variables
gamma = sym('gamma',[l 1],'real');     % the dual variables for the equaltiy
syms t       real;
assume(decision_variables,'real')
if ~isempty(parameters)
    assume(parameters,'real')
end

X_ = [decision_variables; gamma];
KKT_matrix = zeros(n+l,n+l);
KKT_vector = zeros(n+l,1);

%% * First and (second) derivates*
grad_L =0;
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
f_0_tatal = 0 + f_0_cost;

if ~isempty(f_0_cost)
    jac_f_0_cost   = jacobian(f_0_cost,X_); % the first derivate of cost, i.e nebla_f_0
    hess_f_0_cost   = hessian(f_0_cost,X_)  ;
    KKT_matrix = KKT_matrix + hess_f_0_cost ;
    KKT_vector = KKT_vector + -jac_f_0_cost' ;
    grad_L = grad_L + jacobian(f_0_cost,decision_variables)';

end

%2 error terms (cost) (least squares)
for i = 1 : max(size(information_matrix))
    jac_error_cost = jacobian(error{1,i},X_);
    KKT_matrix = KKT_matrix + jac_error_cost'*information_matrix{1,i}*jac_error_cost;
    KKT_vector = KKT_vector + -jac_error_cost'*information_matrix{1,i}*error{1,i};
    f_0_tatal = f_0_tatal + error{1,i}'*information_matrix{1,i}*error{1,i};
    grad_L = grad_L + jacobian(error{1,i},decision_variables)'*information_matrix{1,i}*error{1,i};
end

%--------------------------------------------------------------------------------------

% inequality constraints
if q>0
    f_i = lhs(f_i) - rhs(f_i) ;  % the inequality constraints in form f<=0
end
for i= 1 : q
    error_ineq = f_i(i);
    jac_error_ineq= jacobian(error_ineq,X_);
    slack_inv_i = 1/slack(i);
    Weighted_error = -(lambda(i) + lambda(i) * slack_inv_i * error_ineq + slack_inv_i * 1/t);%lambda_i came from primal dual
    omega_ineq =  lambda(i) * slack_inv_i;

    KKT_matrix = KKT_matrix + jac_error_ineq' * omega_ineq * jac_error_ineq;
    KKT_vector = KKT_vector + jac_error_ineq' * Weighted_error;
    grad_L =  grad_L + (lambda(i) * slack_inv_i * error_ineq + slack_inv_i * 1/t)*jacobian(f_i(i),decision_variables)';

end



% equality cosntraints
omega_eq = [0 1; 1 0];
if l>0
    equality =  lhs(equality) - rhs(equality) ;  % the equality constraints in form g==0
end

for i= 1 : l
    error_eq = [equality(i), gamma(i)]';
    jac_error_eq = jacobian(error_eq,X_);
    KKT_matrix = KKT_matrix + jac_error_eq'*omega_eq*jac_error_eq;
    KKT_vector = KKT_vector + -jac_error_eq'*omega_eq*error_eq;
    grad_L =  grad_L + gamma(i)*jacobian(equality(i),decision_variables)';
end



%% the Newten step


Jacob_f_i       = jacobian(f_i,decision_variables) ;
Jacob_equality = jacobian(equality,decision_variables) ; % the first derivate of equaility, i.e A matrix


Delta_x      = sym('Delta_x',[n 1],'real');  % the update values of the decision variables
if q > 0
Delta_lambda= diag(slack.^-1)*diag(lambda)*(f_i+ Jacob_f_i*Delta_x) + 1/t*diag(slack.^-1);
Delta_slack = -(slack + f_i+ Jacob_f_i*Delta_x) ;
else 
Delta_lambda=[];
Delta_slack = [] ;

end 





%2)###### r_cent
if   ~isempty(f_i)
    r_cent   = diag(lambda)*slack-(1/t)       ;
else
    r_cent = [];
end

%3)###### r_h
if   ~isempty(equality)
    r_pri = equality    ;
else
    r_pri =[];
end

%3)###### r_g
if   ~isempty(f_i)
    r_g = f_i+slack;
else
    r_g =[];
end


%r_t    = [r_dual;r_cent;r_pri;r_g];







%% Preparation for prime daul interior point solver
input= [decision_variables;parameters;slack;lambda;gamma;t];
KKT_matrix_func= sym2func(KKT_matrix,'Vars',input,option) ;
KKT_vector_func= sym2func(KKT_vector,'Vars',input,option);


% residuals vectors
 r_t_func    = sym2func(KKT_vector,'Vars',input,option);
r_dual_func = sym2func(KKT_vector(1:n),'Vars',input,option);
r_pri_func   = sym2func([r_pri;r_g],'Vars',input,option);
% The handle function for cost funtion equality and dynamic cost function
f_func    = sym2func(f_0_tatal,'Vars',input,option); % for evaluation after finishing the calculations

inputForDelta_lambda= [input;Delta_x];
Delta_lambda_func= sym2func(Delta_lambda,'Vars',inputForDelta_lambda,option);
inputForDelta_slack= [input;Delta_x]; % because we use simplified equations 
Delta_slack_func= sym2func(Delta_slack,'Vars',inputForDelta_slack,option);


if   ~isempty(f_i)
    eta_hat = (slack)'*(lambda);
else
    eta_hat =1e-100 ;
end
eta_hat_func= sym2func(eta_hat,'Vars',input,option);
g_func    = sym2func(f_i,'Vars',input,option); % for evaluation after finishing the calculations
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
KKT.f                          = f_0           ;% cost function (symbolic varibales form)
KKT.f_i                        = f_i           ;% inequality function without conversion (symbolic varibales form)
KKT.equality                   = equality           ;% equality function without conversion (symbolic varibales form)
KKT.f_func                     = f_func      ;
KKT.g_func                     = g_func      ; %inequality function of slack variables (symbolic varibales form)
KKT.option                     = option      ; %inequality function of slack variables (symbolic varibales form)
KKT.decision_variables         = decision_variables;
KKT.parameters                 = parameters;
KKT.elimination_flag           = elimination_flag   ;
KKT.algorithm                  ='primal_dual';
KKT.facts                      = [n q l]; %number of the decsion variable;inqualities and equalities
end