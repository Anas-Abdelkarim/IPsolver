function option = options_check(option,decision_variables)

    %options_check(option)
    % option: option are extra options for the solver.
    % details of the option 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 9-a)option.data_recording: takes values 0 and 1 (default) . 1 means 
    % the data record of all variables including the dual variables are
    % saved. This option increase the running time of the code . In case
    % you need to find the optimal solution in shortest time 
    % set option.data_record=0    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 9-b) option.elimination: takes values 0(default) and 1. 1 means
    % the KKT system is based on the reduced KKT matrix (this takes less
    % time). Howeverer, 0 is more robust in numerical calculation. that is
    % why is the default. In case you need to find the optimal solution in 
    % shortest time et option.elimination=0     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % 9-c) option.decision_variables_limits: takes empty (default) or matrix (n,2) n=
    % number of decision varibles. The first column present the lower limit
    % of the decision variables and the second one present the upper limit.
    % the lower limit must be less or equal the upper limit. Otherwise
    % error will appear. The limits are due to the limitation in the
    % functions domain (explicit constraints). For example log(x) x must always larger than zero 
    % we can put this information in this option as option.decision_variables_limits=[0,inf].
    % e.g. 2: consider OP with 4 decision variables this option looks like this:
    % option.decision_variables_limits=[ -10, 10 ;-inf, inf ; -inf, 10; 100, inf ]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 9-d) option.function_structure: '@(in)[matrix]'(default), '@(in){cell}',
    %'@(in1,in2)[matrix]', '@(in1,in2){cell}'. This how we convert the
    %  symbolic variables to handle function. use @(in){cell} or  '@(in1,in2){cell}'
    % if Matlab provides "out of memory" error 
    % 9-10) option.time_save: takes 0 (defualt) and 1; 1 means the
    % solver runs only the core calculation with extra time in recording the data 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 9-e) option.update_dual : takes values 0 and 1(default). This option
    % valid for 'slack-barrier' algorithm. 0 is recommended to be used when
    % all inequalities are linear.       
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 9-f) option.find_feasible_point: takes 0 and 1(default). This option is valid for 
    %'primal-dual-standard' and 'barrier' algorithm. 1 means  
    % to solve an optimization problem to find a start point if the initial 
    % point is not feasible or not given.
 
   
    
    if  isfield(option, 'data_recording')
      switch option.data_recording
         case 0 
             option.data_recording =0;
             disp('option.data_recording=0')
         case 1 
             option.time_save =1 ;
             disp('option.data_recording=1')
         otherwise
             option.time_save =1 ;
             warning('option.data_recording is not correct. The default value is chosen') 
             disp('optiondata_recording=1')
     end
 else
       option.data_recording =1 ;
             disp('option.data_recording=1')
    end   
     
    
      if  isfield(option, 'elimination')
     switch option.elimination
         case 0 
             option.elimination =0;
             disp('option.elimination=0')
         case 1 
             option.elimination =1 ;
             disp('option.elimination=0')
         otherwise
             option.elimination =0 ;
             warning('option.elimination is not correct. The default value is chosen') 
             disp('option.elimination=0')
     end
 else
       option.elimination =0 ;
             disp('option.elimination=0')
      end

  
     if  isfield(option, 'decision_variables_limits')
       if ~isempty(option.decision_variables_limits)
       [r,c] =size(option.decision_variables_limits);
        LL= option.decision_variables_limits(:,1); % lower limit of decion variables
        UL= option.decision_variables_limits(:,2); % upper limit of decion variables        
        check_DVL = UL < LL;
       if r~= length(decision_variables) | c~=2 | sum(check_DVL)~=0
           error('decision_variables_limits are not correct') 
       end 
       inL= ( UL+LL)/2;                                  ; %inside the boundary
       inL(find(isnan(inL)))=1;                            ; % here means no limits
       inL(find(inL==-inf))= UL(find(inL==-inf))-1 ;     ; % means has upper limit
       inL(find(inL== inf))= LL(find(inL== inf))+1 ;     ; % means has lower limit
       option.decision_variables_in =inL; 
       option.Index_decision_variables_up =find(UL<inf) ;
       option.Index_decision_variables_low =find(LL>-inf) ;
       
       end
   else
       option.decision_variables_limits =[] ;
       option.Index_decision_variables_up =[] ;
       option.Index_decision_variables_low =[] ;
       option.decision_variables_in =[] ;
      option.Index_decision_variables_up=[];
       option.Index_decision_variables_low =[]; 
      
   end
   
  
   
    if  isfield(option, 'function_structure')
     switch option.function_structure
         case  '@(in)[matrix]'
             option.function_structure ='@(in)[matrix]';
             disp("option.function_structure='@(in)[matrix]'")
         case '@(in1,in2)[matrix]' 
             option.function_structure ='@(in1,in2)[matrix]' ;
             disp("option.function_structure='@(in1,in2)[matrix]'")
         case '@(in){cell}' 
             option.function_structure ='@(in){cell}';
             disp("option.function_structure='@(in){cell}'")
         case '@(in1,in2){cell}'
             option.function_structure ='@(in1,in2){cell}' ;
             disp("option.function_structure='@(in1,in2){cell}'")
         otherwise
              warning('option.function_structure is not correct. The default value is chosen')
             option.function_structure ='@(in1,in2)[matrix]'  ;
             disp("option.function_structure='@(in1,in2)[matrix]'")
         end
   else
           option.function_structure ='@(in1,in2)[matrix]'  ;
           disp("option.function_structure='@(in1,in2)[matrix]'")
    end
      
   
   
 
 
   if  isfield(option, 'update_dual')
     switch option.update_dual
         case 0 
             option.update_dual =0;
             disp('option.update_dual=0')
         case 1 
             option.update_dual =1 ;
             disp('option.update_dual=1')
         otherwise
             option.update_dual =1 ;
             warning('option.update_dual is not correct. The default value is chosen')
             disp('option.update_dual=1')
     end
 else
       option.update_dual =1 ;
             disp('option.update_dual=1')
   end
   
 
    if  isfield(option, 'find_feasible_point')
     switch option.find_feasible_point
         case 0 
             option.find_feasible_point =0;
             disp('option.find_feasible_point=0')
         case 1 
             option.find_feasible_point =1 ;
             disp('option.find_feasible_point=1')
         otherwise
             option.find_feasible_point =1 ;
             warning('option.find_feasible_point is not correct. The default value is chosen')
             disp('option.find_feasible_point=1')
          end
 else
       option.find_feasible_point =1 ;
             disp('option.find_feasible_point=1')
    end
  
 
    
    option.check_done = 1; % flag to know that this function has been called once.
   
    
end

 
    
    