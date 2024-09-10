function [AL_LS_solver] = solver_AL_LS(KKT,x_start,parameters_subs,warm_start)
tic
% solver_sb(KKT,x_start,parameters_subs,warm_start)
KKT_matrix_func         = KKT.KKT_matrix_func          ;
KKT_vector_func         = KKT.KKT_vector_func          ;
f_0_func                = KKT.f_0_func                 ;
grad_L_func             = KKT.grad_L_func              ;
f_i_func                = KKT.f_i_func                 ;
equality_func           = KKT.equality_func            ;
decision_variables      = KKT.decision_variables       ;
parameters              = KKT.parameters               ;
option                  = KKT.option                   ;
input                   = KKT.input                    ;
n                       = KKT.facts(1)                 ;
q                       = KKT.facts(2)                 ;
l                       = KKT.facts(3)                 ;
function_structure      = option.function_structure    ;
settings                = solver_settings              ;




theta                   = settings.theta               ;
nu                      = settings.nu                  ;
alpha                   = settings.alpha               ;
beta                    = settings.beta                ;
epsilon                 = settings.epsilon             ;
epsilon_feas            = settings.epsilon_feas        ;
iterations_max          = settings.iterations_max      ;

% recodres
if KKT.option.data_recording==1;
    num_rec                    =settings.records_num;
    Newton_calculation_time    =cell(1,num_rec);
    t_record                   =cell(1,num_rec);
    s_record                   =cell(1,num_rec);
    x_record                   =cell(1,num_rec);
    slack_record               =cell(1,num_rec);
    lambda_record              =cell(1,num_rec);
    gamma_record               =cell(1,num_rec);
    Newton_step_record         =cell(1,num_rec);
end
%% intialization
% intialization if x_intitial are not given

if isempty(x_start)
    x_start=settings.x_0*ones(n,1);
end



if isfield (q,'gamma_optimal')
    gamma = warm_start.gamma_optimal;
else
    gamma    = 0*settings.gamma_0*ones(l,1);
end

slack_start= [];
if q>0
    if isfield (warm_start,'slack_optimal')
        slack_start = warm_start.slack_optimal;
    else
        slack_start  = 0*settings.slack_0*ones(q,1);
        %      slack_start = callfunc(g_func,[x_start;parameters_subs],function_structure);
        %      slack_start = settings.slack_0*slack_start;
        %          index_slack_neg= find(slack_start<=.1);
        %       if ~isempty(index_slack_neg)
        %          slack_start(index_slack_neg)= 10;
        %       end
    end
end


if isfield (warm_start,'lambda_optimal')
    lambda     = warm_start.lambda_optimal;
else
    lambda     = 0*settings.lambda_0*ones(q,1);
end

ro_ineq   = ones(q,1) ;
ro_eq     = ones(l,1) ;


lambda = 0*lambda -100;
gamma  = 0*gamma ;
ro_eq   = ro_eq + .2;
ro_ineq   = ro_ineq + .2;
num_iteration =100



% prepeartion of input vector
num_iteration       = 0        ;
x                   = x_start;


f_i    =    callfunc(f_i_func,[x; parameters] ,function_structure);
slack  = max(0,-[lambda./(2*ro_ineq)+f_i; parameters]);


input= [x; parameters; slack; lambda; gamma; ro_ineq; ro_eq];


%% The algorithm
while true
    num_iteration = num_iteration +1 ;
    c_backTracking          = 0 ; % counter of line search (back tracking)
    s = 1;


    f_i    =    callfunc(f_i_func,[x; parameters] ,function_structure);
    slack  = max(0,-[lambda./(2*ro_ineq)+f_i; parameters]);


    input= [x; parameters; slack; lambda; gamma; ro_ineq; ro_eq];

    if KKT.option.data_recording==1;
        %2- ####### data_recording
        x_record{num_iteration}                  = x      ;
        slack_record{num_iteration}              = slack ;
        lambda_record{num_iteration}             = lambda     ;
        gamma_record{num_iteration}              = gamma    ;
        ro_ineq_record{num_iteration}            = ro_ineq    ;
        ro_eq_record{num_iteration}              = ro_eq    ;
        disp(['Augmented Lagrange. Iteration number', num2str(num_iteration)])
    end

    %1-###### Newton step    ##########
    % using unfactorized KKT matrix

    KKT_matrix= callfunc(KKT_matrix_func,input,function_structure);
    KKT_vector= callfunc(KKT_vector_func,input,function_structure);


    Newton_step       = KKT_matrix\KKT_vector ;



    %{
    %2-###### Line search (back-traking) ##########
    s= 1;
    % 4-1 ######## keep slack  positive after updating ###########
    index2= find(lambda+s*Delta_lambda<=0);
    if ~isempty(index2)
        s_safe =  theta*min(-lambda(index2)./Delta_lambda(index2)); % the mininum s that guarantee lambda remain non negative
        s = (min(s,s_safe));
    end

    index_pos_slack= find(slack+s*Delta_slack<=0);
    if ~isempty(index_pos_slack)
        s_safe =  theta*min(-slack(index_pos_slack)./Delta_slack(index_pos_slack)); % the mininum s that guarantee lambda remain non negative
        s = (min(s,s_safe));
        Delta_slack(index_pos_slack)= -theta*slack(index_pos_slack);
    end

    %        % correction of Lambda blocking the movement
    %         if s<s_efficient
    %            if ~isempty(index2)
    %         lambda(index2) = lambda(index2)+ lamda_correction_factor  ;
    %        input= [x;parameters_subs;slack;lambda;gamma;t];
    %        %disp ('lamba small')
    %            end
    %         end

    % 4-2 ######## Implicit constraints backtrack ###########
    if ~isempty(index_upper_limit)||~isempty(index_lower_limit)
        Delta_x= x_in_domain(x,Delta_x,upper_limit,lower_limit,index_upper_limit,index_lower_limit,theta);
    end

    % 4-3 ########### residual backtrack ##############
    r_t = callfunc(r_t_func,input,function_structure);
    while true


        c_backTracking =c_backTracking+1 ;

        x_update     = x       +  s*Delta_x        ; % only proposed update
        lambda_update= lambda  +  s*Delta_lambda   ;
        gamma_update = gamma   +  s*Delta_gamma    ;
        slack_update = slack   +  s*Delta_slack    ;

        input_update =[x_update;parameters_subs;slack_update;lambda_update;gamma_update;t];

        r_t_update = callfunc(r_t_func,input_update,function_structure);
        if norm(r_t_update)>(1-beta*s)*norm(r_t) & s>1e-20
            s=alpha/c_backTracking*s;
        else
            break
        end
    end
    %}



    %4- ###### variables update ##########
    x            = x       + s*Newton_step;

    %5- ###### update other variables  ##########
    f_i    =    callfunc(f_i_func,[x; parameters] ,function_structure);
    equality =   callfunc(equality_func,[x; parameters] ,function_structure);

    lambda       = max(0, lambda  + ro_ineq.*f_i)        ;
    %slack        = max(0,-[lambda./(2*ro_ineq)+f_i; parameters]);

    gamma        = gamma   + 2*ro_eq.*equality            ;
    gamma        = gamma   + 2*equality            ;
    lambda       = max(0, lambda  + 2*f_i)        ;
    slack        = max(0,-[lambda./(2*2)+f_i; parameters]);



    input= [x; parameters; slack; lambda; gamma; ro_ineq; ro_eq];



    %1- ##### termination condition ##########
    grad_L =   callfunc(grad_L_func,input ,function_structure);


    stop = norm(grad_L) <10*epsilon_feas && ...
        norm(equality)  <1*epsilon_feas&&...
        sum(max(0,f_i) >epsilon_feas)==0;
    if stop
        break
    end

    if iterations_max==num_iteration
        warning('The solver has reached the maximum number of iterations')
        break
    end

end

cost_value = callfunc(f_0_func,input,function_structure)                           ;

%% return
AL_LS_solver.x_optimal                  = x                                               ;
AL_LS_solver.cost_value                 = cost_value                                      ;
AL_LS_solver.lambda_optimal             = lambda                                              ;
AL_LS_solver.num_iteration              = num_iteration                                   ;
AL_LS_solver.slack_optimal              = slack                                           ;
AL_LS_solver.gamma_optimal              = gamma                                           ;
if KKT.option.data_recording==1;
    AL_LS_solver.Newton_calculation_time    = array2table([Newton_calculation_time{:}])       ;
    AL_LS_solver.Newton_step_record         = array2table([Newton_step_record{:}])            ;
    AL_LS_solver.x_record                   = array2table([x_record{:}])                      ;
    AL_LS_solver.lambda_record              = array2table([lambda_record{:}])                 ;
    AL_LS_solver.gamma_record               = array2table([gamma_record{:}])                  ;
    AL_LS_solver.slack_record               = array2table([slack_record{:}])                  ;
    AL_LS_solver.search_x_start            = []                                               ;
    % clearvars -except AL_LS_solver
end

AL_LS_solver.solver_time                = toc;

end