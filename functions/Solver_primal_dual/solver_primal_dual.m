function [sb_solver] = solver_primal_dual(KKT,x_start,parameters_subs,warm_start) 
tic
% solver_sb(KKT,x_start,parameters_subs,warm_start) 

KKT_matrix_func         = KKT.KKT_matrix_func          ;
KKT_vector_func         = KKT.KKT_vector_func          ;
Delta_lambda_func       = KKT.Delta_lambda_func        ;
eta_hat_func            = KKT.eta_hat_func             ;
Delta_slack_func        = KKT.Delta_slack_func         ;
r_t_func                = KKT.r_t_func                 ;
r_dual_func             = KKT.r_dual_func              ;
r_pri_func              = KKT.r_pri_func               ;
f_func                  = KKT.f_func                   ;
g_func                  = KKT.g_func                   ;
elimination_flag        = KKT.elimination_flag         ;
option                  = KKT.option                   ;
n                       = KKT.facts(1)                 ;
q                       = KKT.facts(2)                 ;
l                       = KKT.facts(3)                 ;
function_structure      = option.function_structure    ;
settings                = solver_settings              ;
theta                   = settings.theta               ;
tau                     = settings.tau                 ;
nu                      = settings.nu                  ;  
alpha                   = settings.alpha               ;
beta                    = settings.beta                ;
epsilon                 = settings.epsilon             ;
epsilon_feas            = settings.epsilon_feas        ;
iterations_max          = settings.iterations_max      ;
s_efficient             = settings.s_efficient         ;
lambda_correction       = settings.lambda_correction   ;

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
   if isfield (warm_start,'x_optimal')
     x_start = warm_start.x_optimal;
   end
   
   if isempty(x_start)
      if option.find_feasible_point==1;
         x_start=settings.x_0*ones(n,1);
      else 
      disp('ERROR: x_initial is empty. No results are found')
      return
      end
   end
  
 
  %%% check if x_start lambdast satisfies the implicit constraints
   if ~isempty(option.Index_decision_variables_up)
       upper_limit=option.decision_variables_limits(:,2);
       index_upper_limit=  option.Index_decision_variables_up;
       index_check_up= find(x_start(index_upper_limit)>upper_limit(index_upper_limit));
       if sum(index_check_up~=0)
           disp('some variables exceed the upper limit')
           if option.find_feasible_point==1
              x_start(index_upper_limit(index_check_up)) =upper_limit(index_upper_limit(index_check_up));
           end
       end
   end
       
   if ~isempty(option.Index_decision_variables_low)
       lower_limit=option.decision_variables_limits(:,1);
       index_lower_limit=  option.Index_decision_variables_low;
       index_check_low= find(x_start(index_lower_limit)<lower_limit(index_lower_limit));
       if sum(index_check_low~=0)
           disp('some variables exceed the lower limit')
           if option.find_feasible_point==1
              x_start(index_lower_limit(index_check_low)) =lower_limit(index_lower_limit(index_check_low));
           end
       end         
   end
   
   
   if isfield (q,'gamma_optimal')
     gamma = warm_start.gamma_optimal;
   else
     gamma    = settings.gamma_0*ones(l,1);
   end
   
   slack_start= [];
   if q>0 
     if isfield (warm_start,'slack_optimal')
     slack_start = warm_start.slack_optimal;
     else
     %slack_start  = settings.slack_0*ones(q,1);
     slack_start = callfunc(g_func,[x_start;parameters_subs],function_structure); 
     slack_start = settings.slack_0*slack_start;  
         index_slack_neg= find(slack_start<=.1);         
      if ~isempty(index_slack_neg)
         slack_start(index_slack_neg)= 10;         
      end       
     end
   end 

             
     if isfield (warm_start,'lambda_optimal')
     lambda     = warm_start.lambda_optimal;
     else
     lambda     = settings.lambda_0*ones(q,1);
     end     
   
       
% prepeartion of input vector
s                   = 1        ;
num_iteration       = 0        ;
t                   = settings.t_0 ;   
x                   = x_start;
slack               =slack_start;
s=[];

input= [x;parameters_subs;slack;lambda;gamma;t];
%% The algorithm
while true
    num_iteration = num_iteration +1 ;
    c_backTracking          = 0 ; % counter of line search (back tracking)
    c_inTrack     = 0 ; % counter for to keep inside of the interior
 
    
       %1-###### calculate t  ##########
        eta_hat                   = callfunc(eta_hat_func,input,function_structure)     ;
        t                         = nu*q/eta_hat ;               
        input(end)                = t            ;
         
    if KKT.option.data_recording==1; 
     %2- ####### data_recording
    x_record{num_iteration}                  = x      ;
    slack_record{num_iteration}              = slack ;
    lambda_record{num_iteration}                 = lambda     ;
    gamma_record{num_iteration}              = gamma    ;
    s_record{num_iteration}                  = s      ;
     t_record{num_iteration}  = t            ;
    disp(['primal-dual algorithm. Iteration number ', num2str(num_iteration)])
    end
         
    %1-###### Newton step    ##########
    % using unfactorized KKT matrix
     
    KKT_matrix= callfunc(KKT_matrix_func,input,function_structure);
    KKT_vector= callfunc(KKT_vector_func,input,function_structure);

    
    Newton_step       = KKT_matrix\KKT_vector ;

    
    % extract Delta_x, Delta_slack,Delta_lambda and Delta_gamma
       Delta_x     = Newton_step(1    : n ) ;
    if elimination_flag == 0;
     Delta_lambda= Newton_step(n+1  :n+q) ;
     Delta_gamma   = Newton_step(n+q+1:n+q+l) ;   
     Delta_slack   = Newton_step(n+q+l+1:end) ; 
    else % means elimination_flag = 1    
      Delta_gamma= Newton_step(end-l+1:end) ;    
      input_Delta_lambda= [input;  Delta_x];    
      Delta_lambda    = callfunc(Delta_lambda_func,input_Delta_lambda,function_structure) ;        
      input_Delta_slack= [input;  Delta_lambda];
      Delta_slack=callfunc(Delta_slack_func,input_Delta_slack,function_structure);     
   end
    
    
    
      if KKT.option.data_recording==1;  
          if  elimination_flag ==0
            Newton_step_record{num_iteration}   = Newton_step;
          else            
            Newton_step_record{num_iteration}   = [Newton_step;Delta_lambda;Delta_gamma;Delta_slack];         
         end
      end
      
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
      if ~isempty(option.Index_decision_variables_up)
          s= variables_backtracking(x,Delta_x,upper_limit,lower_limit,index_upper_limit,index_lower_limit);     
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
    
     
    
     
    %4- ###### variables update ##########
     x            = x       + s*Delta_x            ;  
     lambda       = lambda      + s*Delta_lambda   ;
     gamma        = gamma  + s*Delta_gamma         ;
     slack        = slack  + s*Delta_slack         ;
      
 
    input        = [x;parameters_subs;slack;lambda;gamma;t];     
    
    %1- ##### termination condition ##########
    r_dual      = callfunc(r_dual_func,input,function_structure)      ;   
    r_pri       = callfunc(r_pri_func,input,function_structure)       ;
    stop = norm(Delta_x) <epsilon_feas && ...
           norm(1/t)<epsilon &&...     
           norm(r_pri)  <epsilon_feas&&...
           norm(r_dual) <sqrt(epsilon_feas);
    if stop
        break
       if KKT.option.data_recording  ==1 ; 
           t_record{num_iteration+1} = t  ;       
           s_record{num_iteration}   = s  ;
       end   
    end
    
   if iterations_max==num_iteration
        break 
        warning('The solver has reached the maximum number of iterations')
    end
   
end
    
cost_value = callfunc(f_func,input,function_structure)                           ;
%% return
sb_solver.x_optimal                  = x                                               ;
sb_solver.cost_value                 = cost_value                                      ;
sb_solver.lambda_optimal             = lambda                                              ;
sb_solver.num_iteration              = num_iteration                                   ; 
sb_solver.slack_optimal              = slack                                           ;
sb_solver.gamma_optimal              = gamma                                           ;
if KKT.option.data_recording==1; 
sb_solver.s_record                   = array2table([s_record{:}])                      ;
sb_solver.Newton_calculation_time    = array2table([Newton_calculation_time{:}])       ;            
sb_solver.Newton_step_record         = array2table([Newton_step_record{:}])            ;
sb_solver.t_record                   = array2table([t_record{:}])                      ; 
sb_solver.x_record                   = array2table([x_record{:}])                      ;
sb_solver.lambda_record              = array2table([lambda_record{:}])                 ;
sb_solver.gamma_record               = array2table([gamma_record{:}])                  ;
sb_solver.slack_record               = array2table([slack_record{:}])                  ;
sb_solver.search_x_start            = []                                               ; 
% clearvars -except sb_solver
end

sb_solver.solver_time                = toc;

end