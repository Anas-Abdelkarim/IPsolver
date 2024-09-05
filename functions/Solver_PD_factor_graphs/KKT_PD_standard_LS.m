function [KKT] = KKT_PD_standard_LS(decision_variables,f_0,f_i,equality,parameters,option);
%KKT_LS(decision_variables,f,g,h,parameters,option);
% PD_= primal- dual st; LS = least squares
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
% f_i is cell of  the inquality constraints
% equality is  cell the inquality constraints
% decision_variables : the decision variables of the optimization  problem
% parameters         : the constants of the optimization problem
%%
% Check the options

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

    %% variables decleration
    n =  length(decision_variables) ; % n is the number of the decision variablesq
    q =  length(f_i);  % 	 is the number of the inqaulity constraints
    l =  length(equality); % l is the number of the eqautilty constraints
    lambda = sym('lambda',[q 1],'real');  % the dual variables for the inequaltiy constriants
    gamma  = sym('gamma',[l 1],'real');     % the dual variables for the equaltiy constriants
    syms t       real;



    %% the Newten step
    if isrow(decision_variables)
        decision_variables = decision_variables';
    end
    X_ = [decision_variables;  lambda; gamma];
    KKT_matrix = zeros(n+q+l,n+q+l);
    KKT_vector = zeros(n+q+l,1);


    %% * First and (second) derivates*

    %1- cost function
    if isstruct(f_0)
        if isfield(f_0, 'cost')
            f_0_cost = f_0.cost;
        else
            f_0_cost =[];
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
    end

    %2 error terms (cost) (least squares)
    for i = 1 : max(size(information_matrix))
        jac_error_cost = jacobian(error{1,i},X_);
        KKT_matrix = KKT_matrix + jac_error_cost'*information_matrix{1,i}*jac_error_cost;
        KKT_vector = KKT_vector + -jac_error_cost'*information_matrix{1,i}*error{1,i};
        f_0_tatal = f_0_tatal + error{1,i}'*information_matrix*error{1,i};
    end

    % inequality constraints
    if q>0
        f_i = lhs(f_i) - rhs(f_i) ;  % the inequality constraints in form f<=0
    end
    for i= 1 : q
        error_ineq = [f_i(i), lambda(i)]';
        jac_error_ineq= jacobian(error_ineq,X_);
        Weighted_error = -[lambda(i); lambda(i)*f_i(i) + 1/t];
        omega_ineq = [0   1
            lambda(i)  f_i(i)];

        KKT_matrix = KKT_matrix + jac_error_ineq' * omega_ineq * jac_error_ineq;
        KKT_vector = KKT_vector + jac_error_ineq' * Weighted_error;
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
    end





    %%%%% Finding feasible start point %%%%%%%
    KKT_x_start=[];
    %% this done by preparing the KKT system of the OP
    if option.find_feasible_point
        if ~isempty (f_i)
            syms level_unique real %
            f_0_x_start= level_unique; % level_unique = k in my forumation
            f_i_x_start = f_i<= level_unique;
            decision_variables_x_start= [ level_unique;decision_variables];
            option.find_feasible_point=0; % to avoid infinite loop
            [KKT_x_start] = KKT_PD_standard_LS(decision_variables_x_start,f_0_x_start,f_i_x_start,[],parameters,option);
            option.find_feasible_point=1; % to return it to its correct value
        end
    end

    %% Preparation for the solver

    input= [decision_variables;parameters;lambda;gamma;t];
    KKT_matrix_func= sym2func(KKT_matrix,'Vars',input,option) ;
    KKT_vector_func= sym2func(KKT_vector,'Vars',input,option);
   
    % surrogate duality gap
    if   ~isempty(f_i)
        eta_hat = (-f_i)'*lambda;
    else
        eta_hat =1e-100 ;
    end
    eta_hat_func= sym2func(eta_hat,'Vars',input,option);
    % residuals vectors
    r_t_func    = sym2func(KKT_vector,'Vars',input,option);
    r_dual_func = sym2func(KKT_vector(1:n),'Vars',input,option);
    r_pri_func  = sym2func(equality,'Vars',input,option);
    % The handle function for cost funtion equality and dynamic cost function
    f_0_func    = sym2func(f_0_tatal,'Vars',input,option); % for evaluation after finishing the calculations
    % grad_f_0_func=sym2func(grad_f_0,'Vars',input,option); % for evaluation after finishing the calculations
    %  equality_func=sym2func(equality,'Vars',input,option);
    if   ~isempty(f_i)
        f_i_func    = sym2func(f_i,'Vars',input,option); % for evaluation after finishing the calculations
    else
        f_i_func=[];
    end

    %% return
    KKT.KKT_matrix_func            = KKT_matrix_func;
    KKT.KKT_vector_func            = KKT_vector_func;
    KKT.eta_hat_func               = eta_hat_func;
    KKT.r_t_func                   = r_t_func    ;
    KKT.r_dual_func                = r_dual_func ;
    KKT.r_pri_func                 = r_pri_func  ;
    KKT.f_0                        = f_0         ;% symbolic varibales for the  cost function
    KKT.f_i                        = f_i     ;% symbolic varibales for the  inequality function without conversion
    KKT.equality                   = equality ;% symbolic variables for the equality==0
    KKT.f_0_func                   = f_0_func    ;
    KKT.f_i_func                   = f_i_func    ;% function of the inquality (with /without slack variables)
    KKT.decision_variables         = decision_variables;
    KKT.parameters                 = parameters        ;
    KKT.elimination_flag           = elimination_flag  ;
    KKT.option                     = option;
    KKT.KKT_x_start                = KKT_x_start;
    KKT.algorithm                  = 'PD_standard_LS';
    KKT.facts                      = [n q l]; %number of the decsion variable;inqualities and equalities
end