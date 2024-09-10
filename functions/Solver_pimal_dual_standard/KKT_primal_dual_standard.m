function [KKT] = KKT_primal_dual_standard(decision_variables,f_0,f_i,equality,parameters,option);
%KKT_sb(decision_variables,f,g,h,parameters,option);.
%This function is run offline. All functions required for the online
%optimization is saved.
% decision_variables : the decision variables of the optimization  problem
% f_0 is the cost function
% f_i is  the inquality constraints
% equality is  the inquality constraints
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
    Delta_x      = sym('Delta_x',[n 1],'real');  % the update values of the decision variables

    %% * First and (second) derivates*
    %cost function
    grad_f_0   = jacobian(f_0,decision_variables)' ; % the first derivate of cost, i.e nebla_f_0
    hess_f_0   = hessian(f_0,decision_variables)  ;
    %inequaility constraints
    if   ~isempty(f_i)
        f = lhs(f_i) - rhs(f_i) ;  % the inequality constraints in form f<=0
        Df       = jacobian(f,decision_variables) ; % the first derivate of inequality, i.e nebla_f^T
        if ~isempty(symvar(Df))  %means some inequality has symbolic variables
            assume(Df,'real')
        end
    else

    end
    %equaility constraints
    if   ~isempty(equality)
        eqaulity             = lhs(equality) - rhs(equality); %the equality constraints in form Ax-b=0
        A                    = jacobian(eqaulity,decision_variables) ; % the first derivate of equality, i.e A matrix

        if ~isempty(symvar(A))  %means some equality has symbolic variables
            assume(A,'real')
        end

    else
        A=[];
        eqaulity=[];
    end
    %% the Newten step
    daulHess_f = zeros(n,n);
    daulGrad_f = zeros(n,n);
    if  ~isempty(f_i)
        for i = 1:q
            daulHess_f =   daulHess_f +lambda(i)*hessian(f(i),decision_variables);
        end



        if elimination_flag==0


            KKT_matrix =[ hess_f_0+daulHess_f   , Df'      , A'
                -diag(lambda)*Df   , -diag(f) , zeros(q,l)
                A          ,zeros(l,q), zeros(l,l) ];


        else

            for i = 1:q
                grad_f{i}  = jacobian(f(i),decision_variables)';
                daulGrad_f =  daulGrad_f -lambda(i)*grad_f{i}*grad_f{i}'/f(i) ;
            end
            H_pd=  hess_f_0 +  daulHess_f +daulGrad_f ;
            KKT_matrix =[ H_pd   ,     A'
                A       , zeros(l,l) ];
        end
    else
        KKT_matrix =[ hess_f_0   ,     A'
            A       , zeros(l,l) ];
    end

    % KKT vector
    %1)###### r_cent

    if   ~isempty(f_i)
        r_cent   = -diag(lambda)*f-(1/t)       ;
    else
        r_cent=[];
    end
    %1)###### r_dual
    r_dual      =  grad_f_0 ;
    if   ~isempty(f_i)
        r_dual     = r_dual+Df'*lambda ;
    end
    if   ~isempty(A)
        r_dual   =  r_dual +A'*gamma ;
    end

    r_pri    =  eqaulity ;

    r_t      = [r_dual; r_cent; r_pri] ;

    Delta_lambda= [];
    if elimination_flag==0 | isempty(f_i)
        KKT_vector= - r_t;
    end
    if elimination_flag==1
        if ~isempty(f_i)
            inequality_inverse=(f).^-1 ;
            KKT_vector= -[r_dual + (Df.*inequality_inverse)'*r_cent; r_pri];
            Delta_lambda= -inequality_inverse.*lambda.*Df*Delta_x+inequality_inverse.*r_cent;
        end
    end
    %%%%% Finding feasible start point %%%%%%%
    KKT_x_start=[];
    %% this done by preparing the KKT system of the OP
    if option.find_feasible_point
        if ~isempty (f_i)
            syms level_unique real %
            f_0_x_start= level_unique;
            f_i_x_start = f<= level_unique;
            decision_variables_x_start= [ level_unique;decision_variables];
            option.find_feasible_point=0; % to avoid infinite loop
            [KKT_x_start] = KKT_primal_dual_standard(decision_variables_x_start,f_0_x_start,f_i_x_start,[],parameters,option);
            option.find_feasible_point=1; % to return it to its correct value
        end
    end

    %% Preparation for the solver
    input= [decision_variables;parameters;lambda;gamma;t];
    KKT_matrix_func= sym2func(KKT_matrix,'Vars',input,option) ;
    KKT_vector_func= sym2func(KKT_vector,'Vars',input,option);
    if   ~isempty(Delta_lambda)
        inputForDelta_lambda= [input;Delta_x];
        Delta_lambda_func= sym2func(Delta_lambda,'Vars',inputForDelta_lambda,option);
    else
        Delta_lambda_func= sym2func(Delta_lambda,'Vars',input,option);
    end
    % surrogate duality gap
    if   ~isempty(f_i)
        eta_hat = (-f)'*lambda;
    else
        eta_hat =1e-100 ;
    end
    eta_hat_func= sym2func(eta_hat,'Vars',input,option);
    % residuals vectors
    r_t_func    = sym2func(r_t,'Vars',input,option);
    r_dual_func = sym2func(r_dual,'Vars',input,option);
    r_pri_func  = sym2func(r_pri,'Vars',input,option);
    % The handle function for cost funtion equality and dynamic cost function
    f_0_func    = sym2func(f_0,'Vars',input,option); % for evaluation after finishing the calculations
    grad_f_0_func=sym2func(grad_f_0,'Vars',input,option); % for evaluation after finishing the calculations
    equality_func=sym2func(equality,'Vars',input,option);
    if   ~isempty(f_i)
        f_i_func    = sym2func(f,'Vars',input,option); % for evaluation after finishing the calculations
    else
        f_i_func=[];
    end

    %% return
    KKT.KKT_matrix_func            = KKT_matrix_func;
    KKT.KKT_vector_func            = KKT_vector_func;
    KKT.Delta_lambda_func          = Delta_lambda_func ;
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
    KKT.algorithm                  = 'primal-dual-standard';
    KKT.facts                      = [n q l]; %number of the decsion variable;inqualities and equalities
end