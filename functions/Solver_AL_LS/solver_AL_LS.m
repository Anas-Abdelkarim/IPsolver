function [AL_LS_solver] = solver_AL_LS(KKT,x_start,parameters_subs,warm_start)
tic
% solver_sb(KKT,x_start,parameters_subs,warm_start)
KKT_matrix_func         = KKT.KKT_matrix_func          ;
KKT_vector_func         = KKT.KKT_vector_func          ;
f_0_func                = KKT.f_0_func                 ;
grad_L_func             = KKT.grad_L_func              ;
f_i_func                = KKT.f_i_func                 ;
equality_func           = KKT.equality_func            ;
parameters              = KKT.parameters               ;
option                  = KKT.option                   ;
n                       = KKT.facts(1)                 ;
q                       = KKT.facts(2)                 ;
l                       = KKT.facts(3)                 ;
slack_flag              = KKT.facts(4)                 ;
function_structure      = option.function_structure    ;
settings                = solver_settings              ;

ro_update               = settings.ro_update           ;
ro_ineq0                = settings.ro_ineq0            ;
ro_eq0                  = settings.ro_eq0              ;
epsilon_feas            = settings.epsilon_feas        ;
iterations_max          = settings.iterations_max      ;
dual_step               = settings.dual_step;
max_num_inner_iterations    = 1;
% recodres
if KKT.option.data_recording==1
    num_rec                    =settings.records_num;
    Newton_calculation_time    =cell(1,num_rec);
    x_record                   =cell(1,num_rec);
    slack_record               =cell(1,num_rec);
    lambda_record              =cell(1,num_rec);
    gamma_record               =cell(1,num_rec);
    Newton_step_record         =cell(1,num_rec);
    ro_ineq_record             =cell(1,num_rec);
    ro_eq_record               =cell(1,num_rec);
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


if isfield (warm_start,'lambda_optimal')
    lambda     = warm_start.lambda_optimal;
else
    lambda     = 0*settings.lambda_0*ones(q,1);
end

ro_eq0 = 2;
ro_ineq0 = 2;
iterations_max  = 120



% prepeartion of input vector
num_outer_iteration       = 0        ;
x                   = x_start;
ro_ineq = ro_ineq0 * ones(q,1) ;
ro_eq   = ro_eq0 * ones(l,1) ;


f_i    =    callfunc(f_i_func,[x; parameters] ,function_structure);
if slack_flag
    active  = double(f_i>-.5*lambda./ro_ineq);
else
    active  = double(f_i>0);
end
input= [x; parameters; active; lambda; gamma; ro_ineq; ro_eq];
num_iterations= 0 ;
%% The algorithm
while true

    num_outer_iteration = num_outer_iteration +1 ;

    if KKT.option.data_recording==1
        %2- ####### data_recording
        x_record{num_outer_iteration}                  = x      ;
        lambda_record{num_outer_iteration}             = lambda    ;
        gamma_record{num_outer_iteration}              = gamma    ;
        ro_ineq_record{num_outer_iteration}            = ro_ineq    ;
        ro_eq_record{num_outer_iteration}              = ro_eq    ;
    end

    %1-###### Newton step    ##########
    % using unfactorized KKT matrix
    Newton_step = 1;


    for k = 1 : max_num_inner_iterations
        if   k> iterations_max || num_iterations >=iterations_max
            break
        end

        num_iterations = num_iterations +1;
        KKT_matrix = callfunc(KKT_matrix_func,input,function_structure);
        KKT_vector = callfunc(KKT_vector_func,input,function_structure);


        Newton_step       = (KKT_matrix+0*eye(n,n))\KKT_vector;



        %4- ###### variables update ##########
        s  =1;    x  =  x  + s*Newton_step;
        lambda
        x

        if KKT.option.data_recording==1
            disp(['Augmented Lagrange. Iteration number', num2str(num_iterations-1)])
        end
        if sum(isnan(Newton_step)>0)
            num_outer_iteration
            iterations_max
        end
        input= [x; parameters; active; lambda; gamma; ro_ineq; ro_eq];

        if norm(Newton_step)< .01
            break
        end
    end
    num_iteration_inner(num_outer_iteration) = k;


    %5- ###### update other variables  ##########
    f_i    =    callfunc(f_i_func,[x; parameters] ,function_structure);
    equality =   callfunc(equality_func,[x; parameters] ,function_structure);

    lambda       = max(0, lambda  + 2*ro_ineq.*f_i)   ;

    if slack_flag
        active  = double(f_i>-.5*lambda./ro_ineq);
    else
        active  = double(f_i>0);
    end
    gamma        = gamma + 2*ro_eq.*equality         ;

    active_record{num_outer_iteration}              = active ;





    input= [x; parameters; active; lambda; gamma; ro_ineq; ro_eq];



    %1- ##### termination condition ##########
    grad_L =   callfunc(grad_L_func,input ,function_structure);

    stop = norm(grad_L) <epsilon_feas && ...
        norm(equality)  <epsilon_feas&&...
        norm(Newton_step) <epsilon_feas && ...
        sum(max(0,f_i) >epsilon_feas)==0;
    if stop
        if KKT.option.data_recording==1
            %2- ####### data_recording
            x_record{num_outer_iteration}                  = x      ;
            gamma_record{num_outer_iteration}              = gamma    ;
            ro_ineq_record{num_outer_iteration}            = ro_ineq    ;
            ro_eq_record{num_outer_iteration}              = ro_eq    ;
        end
        break
    end

    if num_iterations>=iterations_max
        warning('The solver has reached the maximum number of iterations')
        if KKT.option.data_recording==1
            %2- ####### data_recording
            x_record{num_outer_iteration}                  = x      ;
            lambda_record{num_outer_iteration}             = lambda    ;
            gamma_record{num_outer_iteration}              = gamma    ;
            ro_ineq_record{num_outer_iteration}            = ro_ineq    ;
            ro_eq_record{num_outer_iteration}              = ro_eq    ;
            disp(['Augmented Lagrange. Iteration number', num2str(num_outer_iteration)])
        end
        break
    end





end




cost_value = callfunc(f_0_func,input,function_structure)                           ;

%% return
AL_LS_solver.x_optimal                  = x                                               ;
AL_LS_solver.cost_value                 = cost_value                                      ;
AL_LS_solver.lambda_optimal             = lambda                                          ;
AL_LS_solver.num_iteration              = num_iterations                                  ;
AL_LS_solver.slack_optimal              = active                                          ;
AL_LS_solver.gamma_optimal              = gamma                                           ;

if KKT.option.data_recording==1
    AL_LS_solver.Newton_calculation_time    = array2table([Newton_calculation_time{:}])       ;
    AL_LS_solver.Newton_step_record         = array2table([Newton_step_record{:}])            ;
    AL_LS_solver.x_record                   = array2table([x_record{:}])                      ;
    AL_LS_solver.lambda_record              = array2table([lambda_record{:}])                 ;
    AL_LS_solver.gamma_record               = array2table([gamma_record{:}])                  ;
    AL_LS_solver.slack_record               = array2table([active_record{:}])                  ;
    AL_LS_solver.ro_ineq_record             = array2table([ro_ineq_record{:}])                ;
    AL_LS_solver.ro_eq_record               = array2table([ro_eq_record{:}])                  ;
    AL_LS_solver.search_x_start             = []                                              ;
    % clearvars -except AL_LS_solve
end

AL_LS_solver.solver_time                = toc;

end