function [bar_solver] = solver_barrier(KKT,x_start,parameters_subs,warm_start)
%solver_bar(KKT,x_start,parameters_subs,warm_start)
tic
bar_solver = [] ;
KKT_matrix_func       = KKT.KKT_matrix_func         ;
KKT_vector_func       = KKT.KKT_vector_func         ;
f_0_func              = KKT.f_0_func                ;
f_i_func              = KKT.f_i_func                ;
decrement_func        = KKT.decrement_func          ;
grad_L_func           = KKT.grad_L_func             ;
equality_func          =KKT.equality_func           ;
option                = KKT.option                  ;
n                     = KKT.facts(1)                ; %number of the decsion variable
q                     = KKT.facts(2)                ; %number of the inqualities
l                     = KKT.facts(3)                ; %number of the  equalities
settings              = solver_settings             ;
nu                    = settings.nu                 ;
theta                 = settings.theta              ;
alpha                 = settings.alpha              ;
beta                  = settings.beta               ;
epsilon               = settings.epsilon            ;
epsilon_feas          = settings.epsilon_feas        ;
iterations_max        = settings.iterations_max     ;

accept_warm_start_ineq_flag   = 0;
accept_warm_start_eq_flag   = 0;
search_x_start_flag         = 0;
function_structure          = option.function_structure    ;

% recodres
if KKT.option.data_recording==1;
    num_rec                    = 5*settings.records_num;
    t_record                   =cell(1,num_rec);
    s_record                   =cell(1,num_rec);
    x_record                   =cell(1,num_rec);
    Newton_step_record         =cell(1,num_rec);
    mu_record                  =cell(1,num_rec);
    gamma_record               =cell(1,num_rec);
    slack_record               =cell(1,num_rec);
end
%% intialization


% check if warm_start statisfies the inequality
if isfield (warm_start,'x_optimal')
    %% check the feasibility of warm_start if incase the slack variables are not used
    if  ~isempty(KKT.f_i) % we check if there is inequality to satisfy, the constraints on the salck are ignored
        input =[warm_start.x_optimal;parameters_subs];
        if sum(callfunc(f_i_func,input,function_structure)>0)==0
            x_start = warm_start.x_optimal;
            accept_warm_start_ineq_flag= 1;
        else
            warning('Warm start gives infeasible x_start; the warm start is ingonred')
        end
    else % means we accept the warm start because there is slack variables or no inequality constraints
        x_start = warm_start.x_optimal;
        accept_warm_start_ineq_flag= 1;
    end
end

if isempty(x_start) % means we do not have x_start at all
    if    option.find_feasible_point==1
        x_start= ones(n,1);
    else
        warning('x_initial is not given')
        return
    end
end


if ~isempty(x_start) & ~accept_warm_start_ineq_flag
    %% check the feasibility of x_start
    if ~isempty(KKT.f_i) % we check if there is inequality to satisfy, the constraints on the salck are ignored
        input =[x_start;parameters_subs];
        if sum(callfunc(f_i_func,input,function_structure)>=0)~=0
            warning('x_initial is not feasible')
            if option.find_feasible_point==1 % try to find x_strat because
                disp('WARNING: x_initial is empty. We attempt to find a x_initial')
                disp('WARNING: To cancel this option set option.find_feasible_point=0')
                search_x_start_flag=1 ;
                if ~isempty(KKT.KKT_x_start.x_start_search)
                    A= callfunc(KKT.KKT_x_start.x_start_search.A_func,input,function_structure);
                    b= callfunc(KKT.KKT_x_start.x_start_search.b_func,input,function_structure);
                    x_start_search = A\b ;
                else
                    x_start_search= ones(n,1);
                end
                input     = [x_start_search;parameters_subs];
                level_unique=  max(callfunc(f_i_func,input,function_structure)) +10;
                x_start_search = [level_unique; x_start_search];
            end
        end
    end




    %%%%%%%%% call the solver to find x_start %%%%%%%%
    search_x_start=[];
    if search_x_start_flag==1
        KKT.KKT_x_start.option.search_x_start_flag=1;
        [search_x_start] = solver_barrier(KKT.KKT_x_start,x_start_search,parameters_subs,[])
        x_start= search_x_start.x_optimal(2:end);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%% check if x_start must satisfies the implicit constraints
    index_upper_limit=  option.Index_decision_vars_up;
    index_lower_limit=  option.Index_decision_vars_low;

    if ~isempty(index_upper_limit)||~isempty(index_lower_limit)
        lower_limit=option.decision_variables_limits(:,1);
        upper_limit=option.decision_variables_limits(:,2);
    end

    %%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Initialization
    num_iteration   = 0         ;
    num_iteration_Newton=0 ;
    num_inner_iteration = 0      ;
    c_num_inner = 0;
    %% Initialization
    if isfield (warm_start,'gamma_optimal')
        gamma = warm_start.gamma_optimal;
    else
        gamma    = settings.gamma_0*ones(l,1);
    end

    s=[];
    t                = settings.t_0        ;
    x= x_start                             ;
    input        = [x;parameters_subs;gamma;t];


    nu = 30
    t = .1;
    t_array= {t, 2 ,300,1200};


    %% The algorithm
    while true
        num_inner_iteration  = 0         ;
        c_backTracking    = 0 ; % counter of line search
        c_inTrack=0;
        c_num_inner = c_num_inner + 1;
          

        disp(['Barrier algorithm. Iteration number ', num2str(num_iteration)])




        %1-###### Beginning of the inner loop  ##########

        while true
                    num_iteration = num_iteration +1 ;


            num_inner_iteration  = num_inner_iteration+1         ;

            %%%% 1- Calculate the Newton step

            num_iteration_Newton= num_iteration_Newton+1;
            KKT_matrix= callfunc(KKT_matrix_func,input,function_structure);
            KKT_vector= callfunc(KKT_vector_func,input,function_structure);
            Newton_step       = KKT_matrix\KKT_vector ;

            % 1- ####### data_recording
            if KKT.option.data_recording==1 ;
                x_record{num_iteration}                = x      ;
                gamma_record{num_iteration}            = gamma  ;
                s_record{num_inner_iteration,num_iteration} = s              ;
                Newton_step_record{num_iteration_Newton} = Newton_step    ;
            end

            % extact Delta_x,
            Delta_x  = Newton_step(1:n)     ;
            Delta_gamma    = Newton_step(n+1:end)   ;







            %%% 2- stopping criterion for the inner loop
            inputForDecrement= [input; Delta_x];
            decrement=callfunc(decrement_func,inputForDecrement,function_structure);

            % if norm(decrement) <= epsilon_feas
            %  break
            % end
            s= 1;
            %
            % %%%%% 3- Line search
            %

            % % 3-A ######## Implicit constraints backtrack ###########
            % if ~isempty(index_upper_limit)||~isempty(index_lower_limit)
            %   Delta_x= x_in_domain(x,Delta_x,upper_limit,lower_limit,index_upper_limit,index_lower_limit,theta);
            % end
            %
            % %%%%%%% 3-B  check if we are still inside the interior

            if ~isempty(f_i_func)
                while true
                    c_inTrack = c_inTrack+1 ;
                    x_update     = x       + s*Delta_x     ;
                    input_update = [x_update;parameters_subs;gamma;t];
                    f_i_update  =callfunc(f_i_func,input_update,function_structure);
                    index_out_interior{c_inTrack}= find(f_i_update>=-0.000000000001);

                    if ~isempty(index_out_interior{c_inTrack})& s>1e-20
                        s=alpha*s;

                        %disp('comming back to feasibility')
                    else
                        break
                    end
                end
            end

            
            
            


            %4- ###### variables update ##########

            x      = x     + s*Delta_x                    ;
            gamma  = gamma + s*Delta_gamma                ;

            input  = [x;parameters_subs;gamma;t]          ;

            if   isfield(option,'search_x_start_flag')
                if x(1)<0
                    break
                end
            end

            % break to test my updating t
            grad_L=callfunc(grad_L_func,input,function_structure);
            if norm(grad_L) <= epsilon_feas *100
                num_inner_iteration_record{c_num_inner} =   num_inner_iteration
                break

            end
        end


        % ###### end of the inner loop  ##########


        t_record{num_iteration} = t;

        %3- ##### Stopping criterion for barrier method ##########


        if  1/t<=epsilon     ||q==0
            f_i  =callfunc(f_i_func,input,function_structure);
            r_pri=callfunc(equality_func,input,function_structure);
            stop = norm(grad_L) <epsilon_feas && ...
                norm(r_pri)  <epsilon_feas&&...
                sum(max(0,f_i) >epsilon_feas)==0  ;
            if stop

                break
            end
        end
        %    tau=.01;nu=50;   to test my updating t
        %    t                        = t*(s*exp(-tau*norm(Delta_x))*nu+1) ;
        if t<t_array{end}
            t = t_array{c_num_inner+1};
        end
        input(end) = t         ;


    end
    
    cost_value =callfunc(f_0_func,input,function_structure)                          ;
    %% return
    bar_solver.x_optimal                  = x                                           ;
    bar_solver.cost_value                 = cost_value                                  ;
    bar_solver.gamma_optimal              = gamma                                       ;
    bar_solver.num_iteration              = num_iteration                               ;
    bar_solver.num_iteration_Newton       = num_iteration_Newton                        ;
    if KKT.option.data_recording==1                                                     ;
        bar_solver.s_record                   = array2table([s_record{:}])                  ;
        bar_solver.Newton_step_record         = array2table([Newton_step_record{:}])        ;
        bar_solver.t_record                   = array2table([t_record{:}])                  ;
        bar_solver.x_record                   = array2table([x_record{:}])                  ;
        bar_solver.gamma_record               = array2table([gamma_record{:}])              ;
        bar_solver.num_inner_iteration_record = array2table([num_inner_iteration_record{:}]);
        bar_solver.search_x_start             = search_x_start                              ;
    end
    bar_solver.solver_time                = toc                                         ;
end