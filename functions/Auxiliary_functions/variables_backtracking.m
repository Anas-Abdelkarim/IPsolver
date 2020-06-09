function s= variables_backtracking(x,Delta_x,upper_limit,lower_limit,index_upper_limit,index_lower_limit) 
         s=1;
         index_check_up= find(x(index_upper_limit)+Delta_x(index_upper_limit)>upper_limit(index_upper_limit));
         if sum(index_check_up~=0)
           %  disp('some variables exceed the upper limit')
           if option.find_feasible_point==1
               s_safe =  .99*min(x(index_check_up)./Delta_x(index_check_up)); % the mininum s that guarantee lambda remain non negative  
               s      = (min(s,s_safe));
               % Delta_x(index_upper_limit(index_check_up)) =.95*x(index_upper_limit(index_check_up));
           end
        end
      end
     
      
      if ~isempty(option.Index_decision_variables_low)
         index_check_low=find(x(index_lower_limit)+Delta_x(index_lower_limit)<lower_limit(index_lower_limit));
         if sum(index_check_low~=0)
            if option.find_feasible_point==1
                s_safe =  .99*min(-x(index_check_low)./Delta_x(index_check_low)); % the mininum s that guarantee lambda remain non negative  
                s      = (min(s,s_safe));
               % Delta_x(index_lower_limit(index_check_low)) =-.95*x(index_lower_limit(index_check_low));
            end
         end         
 end