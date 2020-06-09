function [pd_solver] = solver_primal_dual_standard(KKT,x_start,parameters_subs,warm_start); 
tic

pd_solver               = []                           ;
KKT_matrix_func         = KKT.KKT_matrix_func          ;
KKT_vector_func         = KKT.KKT_vector_func          ;
Delta_lambda_func       = KKT.Delta_lambda_func     ;
eta_hat_func            = KKT.eta_hat_func             ;
r_t_func                = KKT.r_t_func                 ; 
r_dual_func             = KKT.r_dual_func              ;
r_pri_func              = KKT.r_pri_func               ;    
f_i_func                = KKT.f_i_func                 ;
f_0_func                = KKT.f_0_func                 ;
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
accept_warm_start_flag  = 0                            ;
search_x_start_flag     = 0                            ;
% recodres
if KKT.option.data_recording==1; 
num_rec                    =settings.records_num;
t_record                   =cell(1,num_rec);
s_record                   =cell(1,num_rec);
x_record                   =cell(1,num_rec);
Newton_step_record         =cell(1,num_rec); 
lambda_record              =cell(1,num_rec);
mu_record                  =cell(1,num_rec);
gamma_record               =cell(1,num_rec);
end
    
%% x_start : check of x_start 

 %%%check the feasibility of warm_start if  the slack variables are not used
   if isfield (warm_start,'x_optimal')
       if  ~isempty(KKT.f_i) % we check if there is inequality to satisfy, the constraints on the salck are ignored
         input =[warm_start.x_optimal;parameters_subs];
             if sum(callfunc(f_i_func,input,function_structure)>0)==0
             x_start = warm_start.x_optimal;  
             accept_warm_start_flag= 1; 
             else
             warning('Warm start gives infeasible x_start; the warm start is ingonred')
             end
       else % means we accept the warm start because there is no inequality constraints
          x_start = warm_start.x_optimal;
          accept_warm_start_flag= 1; 
       end
   end
   
   
   
    if isempty(x_start) % means we do not have x_start at all
    if    option.find_feasible_point==1
          x_start=settings.x_0*ones(n,1);
    else
         warning('x_initial is not given')
         return       
    end
    end
      
   
   
   if ~isempty(x_start) & ~accept_warm_start_flag &  ~isfield(option,'search_x_start_flag') 
       %% check the feasibility of x_start
       if  ~isempty(KKT.f_i) % we check if there is inequality to satisfy, the constraints on the salck are ignored
         input =[x_start;parameters_subs];
             if sum(callfunc(f_i_func,input,function_structure)>=0)~=0
                  warning('x_initial is not feasible')
                  if option.find_feasible_point==1 % try to find x_strat because 
                      disp('WARNING: x_initial is empty. We attempt to find a x_initial')  
                      disp('WARNING: To cancel this option set option.find_feasible_point=0')
                      search_x_start_flag=1 ;
                      if ~isempty(option.decision_variables_in)
                      x_start_search =  option.decision_variables_in;  
                      else
                      x_start_search= settings.x_0*ones(n,1);  
                      end
                      input     = [x_start_search;parameters_subs]; 
                      level_unique=  max(callfunc(f_i_func,input,function_structure)) +10;
                      x_start_search = [level_unique; x_start_search];                       
                  end
             end
       end
   end
   
    %%%%%%%% check if there are implicit constraints
   if ~isempty(option.Index_decision_variables_up)
       upper_limit=option.decision_variables_limits(:,2);
       index_upper_limit=  option.Index_decision_variables_up;
   end
       
   if ~isempty(option.Index_decision_variables_low)
       lower_limit=option.decision_variables_limits(:,1);
       index_lower_limit=  option.Index_decision_variables_low;
   end   
  %%%%%%%%%%%%%%%%%     
   
   %%%%%%%%% call the solver to find x_start %%%%%%%%
   search_x_start=[];
   if search_x_start_flag==1
     KKT.KKT_x_start.option.search_x_start_flag=1;
     [search_x_start] = solver_primal_dual_standard(KKT.KKT_x_start,x_start_search,parameters_subs,[]);  
      x_start= search_x_start.x_optimal(2:end);             
   end
   
    %%%%%%%% check there are implicit constraints
       if ~isempty(option.Index_decision_variables_up)
       upper_limit=option.decision_variables_limits(:,2);
       index_upper_limit=  option.Index_decision_variables_up;
       end
       
       if ~isempty(option.Index_decision_variables_low)
       lower_limit=option.decision_variables_limits(:,1);
       index_lower_limit=  option.Index_decision_variables_low;
       end   
    %%%%%%%%%%%%%%%%% 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %% Initialization 
   if isfield (warm_start,'gamma_optimal')
     gamma = warm_start.gamma_optimal;
   else
     gamma    = settings.gamma_0*ones(l,1);
   end
   
   if isfield (warm_start,'lambda_optimal')
     lambda = warm_start.lambda_optimal;
   else
     lambda     = settings.lambda_0*ones(q,1);
   end
   
   
   
      
   num_iteration   = 0           ;
   s           =[];     
   t           = settings.t_0 ;   
   x           = x_start;
  input        = [x_start;parameters_subs;lambda;gamma;t]; 
   
      
%% The algorithm
while true
    num_iteration = num_iteration +1 ;
    c_backTracking          = 0 ; % counter of line search (back tracking)
    c_inTrack     = 0 ; % counter for to keep inside of the interior
   
    %1-###### calculate t  ##########
     eta_hat                   = callfunc(eta_hat_func,input,function_structure)     ;
     t                         = nu*q/eta_hat ;
     input(end)               = t            ;
  
    
    
    
   %2- ####### data_recording
    if KKT.option.data_recording==1; 
    t_record{num_iteration}                = t          ;
    x_record{num_iteration}                = x          ;
    lambda_record{num_iteration}           = lambda     ;
    gamma_record{num_iteration}            = gamma      ;
    s_record{num_iteration}                = s          ;
    end
    disp(['Primal-dual standard algorithm. Iteration number ', num2str(num_iteration)])

  
     
    %2-###### Newton step    ##########
    % using unfactorized KKT matrix
     
    KKT_matrix= callfunc(KKT_matrix_func,input,function_structure);
    KKT_vector= callfunc(KKT_vector_func,input,function_structure);
%    
%     Newton_step_linsolve= linsolve(KKT_matrix,KKT_vector);
%     Newton_step_linsolve_time{num_iteration}=toc;
%    
    
    Newton_step       = KKT_matrix\KKT_vector ;
   
  
    
    % extact Delta_x, Delta_gamma and     Delta_lambda
    Delta_x     = Newton_step(1    : n ) ;
   
    if elimination_flag==0
    Delta_lambda= Newton_step(n+1  :n+q) ;
    Delta_gamma    = Newton_step(end-l+1:end) ;
    else
    Delta_gamma = Newton_step(end-l+1:end) ;
    input_lambda= [input;  Delta_x];
    Delta_lambda=callfunc(Delta_lambda_func,input_lambda,function_structure);
    end
    
    
    
    
    Newton_step_record{num_iteration}        = Newton_step ;
  
     %3-###### Line search (back-traking) ##########
     s= 1;
     % 3-1 ######## keep lambda  positive after updating ###########
      index2= find(lambda+s*Delta_lambda<=0);
        
      if ~isempty(index2)
       s_safe =  theta*min(-lambda(index2)./Delta_lambda(index2)); % the mininum s that guarantee lambda remain non negative  
       s = (min(s,s_safe));
      end
        
     % correction of Lambda blocking the movement
        if s<s_efficient 
           if ~isempty(index2)
        lambda(index2) = lambda(index2)+ lambda_correction  ;
       input        = [x;parameters_subs;lambda;gamma;t]; 
          %disp ('lamba small')
           end  
        end
        
        
      % 3-2 ######## Implicit constraints backtrack ###########
      if ~isempty(option.Index_decision_variables_up)
          s= variables_backtracking(x,Delta_x,upper_limit,lower_limit,index_upper_limit,index_lower_limit);     
      end
   
    % 3-3 ########### residual backtrack ##############
  
    r_t= callfunc(r_t_func,input,function_structure);
     
    while true
    c_backTracking =c_backTracking+1 ;
    
    lambda_update= lambda  + s*Delta_lambda  ;% only proposed update
    gamma_update = gamma      + s*Delta_gamma      ; 
    x_update     = x       + s*Delta_x       ; 
   
    input_update =[x_update;parameters_subs;lambda_update;gamma_update;t];
   
    r_t_update = callfunc(r_t_func,input_update,function_structure);
     if norm(r_t_update)>(1-beta*s)*norm(r_t)  & s>0
        s=alpha/c_backTracking*s;
         else                  
     break       
     end 
    end
     
     
     
%     3-4 ######## check if we are still inside the interior 

if ~isempty(f_i_func)
     while true
     c_inTrack = c_inTrack+1 ;
     
     f_i_update  =callfunc(f_i_func,input_update,function_structure);
     index_out_interior{num_iteration,c_inTrack}= find(f_i_update>=-0.000000000001);
     
      if ~isempty(index_out_interior{num_iteration,c_inTrack})  & s>0
       s=alpha/ c_inTrack*s  ;       
       lambda_update = lambda  + s*Delta_lambda  ;% only proposed update
       gamma_update  = gamma   + s*Delta_gamma      ; 
       x_update      = x       + s*Delta_x       ; 
     
      input_update =[x_update;parameters_subs;lambda_update;gamma_update;t];
       %disp('comming back to feasibility')
      else 
          break
      end
     end
end

     
    %4- ###### variables update ##########
        
    lambda       = lambda  + s*Delta_lambda              ;
    gamma        = gamma      + s*Delta_gamma                   ;
    x             = x       + s*Delta_x             ;      
    input        = [x;parameters_subs;lambda;gamma;t]; 
    
 
   
      %5- ##### termination condition ##########
    r_daul      = callfunc(r_dual_func,input,function_structure)      ;     
    r_pri       = callfunc(r_pri_func,input,function_structure)       ;
    eta_hat     = callfunc(eta_hat_func,input,function_structure)     ;
    stop = norm(r_daul) <epsilon_feas &&...
           norm(r_pri)  <epsilon_feas &&...
           norm(eta_hat)<epsilon             ; 
    if stop
        break
        if KKT.option.data_recording == 1 ; 
           t_record{num_iteration+1} = t  ;       
           s_record{num_iteration}   =s  ;
       end   
    end
  
    if   isfield(option,'search_x_start_flag')  
    if x(1)<0
        break
    end
    end
    
       
   if iterations_max==num_iteration
        break 
        warning('The solver has reached the maximum number of iterations')
    end
   
    
    
  
 
end

    
cost_value = callfunc(f_0_func,input,function_structure)                            ;
%% return
 pd_solver.x_optimal                  = x                         ;
 pd_solver.gamma_optimal              = gamma                          ;
 pd_solver.cost_value                 = cost_value                                     ;
 pd_solver.num_iteration              = num_iteration                                  ; 
if KKT.option.data_recording==1; 
 pd_solver.s_record                   = array2table([s_record{:}])                     ;
 pd_solver.Newton_step_record         = array2table([Newton_step_record{:}])           ;
 pd_solver.t_record                   = array2table([t_record{:}])                     ;
 pd_solver.x_record                   = array2table([x_record{:}])                     ;
 pd_solver.gamma_record               = array2table([gamma_record{:}])                 ; 
 pd_solver.lambda_record              = array2table([lambda_record{:}])                ;
 pd_solver.search_x_start             = search_x_start                                 ; 
end
 pd_solver.lambda_optimal             = lambda                                         ;
 pd_solver.solver_time                = toc;
 %pd_solver = pd_solver (~ cellfun (@isempty, {pd_solver.places})); % delete empty struct 
end