function [KKT] = KKT_barrier_LS(decision_variables,f_0,f_i,equality,parameters,option)
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
    if ~exist('option','var')
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
    
    if option.update_dual
        update_dual = 1;
    else
        update_dual = 0;
    end 


end

%% variables decleration
n =  length(decision_variables) ; % n is the number of the decision variablesq
q =  length(f_i);  % 	 is the number of the inqaulity constraints
l =  length(equality); % l is the number of the eqautilty constraints
gamma  = sym('gamma',[l 1],'real');     % the dual variables for the equaltiy constriants
syms t       real;

assume(decision_variables,'real')
if ~isempty(parameters)
    assume(parameters,'real')
end
%% the Newten step
if isrow(decision_variables)
    decision_variables = decision_variables';
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


 % inequality constraints
    if q>0
        f_i = lhs(f_i) - rhs(f_i) ;  % the inequality constraints in form f<=0
    end
    for i= 1 : q
        error_ineq = f_i(i);
        jac_error_ineq= jacobian(error_ineq,X_);
        Weighted_error = 1/f_i(i)/t;
        omega_ineq = 1/f_i(i)^2/t; % assuming no second order 

        KKT_matrix = KKT_matrix + jac_error_ineq' * omega_ineq * jac_error_ineq;
        KKT_vector = KKT_vector + jac_error_ineq' * Weighted_error;
        grad_L =  grad_L - jacobian(error_ineq,decision_variables)' * Weighted_error;

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
        if update_dual
            KKT_vector = KKT_vector + -jac_error_eq'*omega_eq*error_eq;
        else
            Weighted_error = [0;0];
            KKT_vector = KKT_vector + jac_error_eq' * Weighted_error;
        end
        grad_L =  grad_L + gamma(i)*jacobian(equality(i),decision_variables)';
    end






Delta_x_bar      = sym('Delta_x_bar',[n 1],'real');  % the update values of the decision variables
decrement=  Delta_x_bar'*KKT_vector(1:n) ;



%%%%% Finding feasible start point %%%%%%%
KKT_x_start=[];
%% this done by preparing the KKT system of the OP
if option.find_feasible_point
     if ~isempty (f_i)
        syms level_unique real %
        f_0_x_start = level_unique;
        f_i_x_start = f_i <= level_unique;
        decision_variables_x_start = [level_unique;decision_variables];
        option.find_feasible_point = 0; % to avoid infinite loop
        KKT_x_start = KKT_barrier_LS(decision_variables_x_start,f_0_x_start,f_i_x_start,[],parameters,option);
    end 
    option.find_feasible_point=1; % to return it to its correct value
end




%% Preparation for the interior point solver

input= [decision_variables;parameters;gamma;t];
% the newton step
KKT_matrix_func= sym2func(KKT_matrix,'Vars',input,option) ;
%KKT_matrix_factorized.L= matlabFunction(L,'Vars',input,option);
%KKT_matrix_factorized.U= matlabFunction(U,'Vars',input,option);
%KKT_matrix_factorized.P= P                             ;
KKT_vector_func= sym2func(KKT_vector,'Vars',input,option);
grad_L_func  = matlabFunction(grad_L,'Vars',input); % for evaluation after finishing the calculations


% The cost function
f_0_func    = sym2func(f_0_tatal,'Vars',input,option); % for evaluation after finishing the calculations
%f_0_bar_func      = sym2func(f_0_bar,'Vars',input,option); % the approximated cost function
f_i_func  = sym2func(f_i,'Vars',input,option); % for backtracking
equality_func    = sym2func(equality,'Vars',input,option); %
inputForDecrement= [input;Delta_x_bar];
decrement_func= sym2func(decrement,'Vars',inputForDecrement,option);


%% return
KKT.KKT_matrix_func             = KKT_matrix_func;
%KKT.KKT_matrix_factorized      = KKT_matrix_factorized;
KKT.KKT_vector_func             = KKT_vector_func;
KKT.f_0_func                    = f_0_func   ; % for evaluation after finishing the calculations
KKT.f_i_func                    = f_i_func   ;  % for evaluation after finishing the calculations
KKT.grad_L_func                 = grad_L_func    ;
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
KKT.algorithm                   = 'barrier_LS'      ;
KKT.facts                       = [n q l update_dual ]; %number of the decsion variable;inqualities and equalities
end