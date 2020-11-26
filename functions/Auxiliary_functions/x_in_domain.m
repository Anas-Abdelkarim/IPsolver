function output= x_in_domain(x,Delta_x,upper_limit,lower_limit,index_upper_limit,index_lower_limit,theta,in_domain)
% this function check if the decision variables inside the domain of the
% optimization problem 


 if  isempty(Delta_x) % means we prepare for x_start
     output=x;
     if ~isempty(index_upper_limit)
          index_check_up= find(x(index_upper_limit)>upper_limit(index_upper_limit));
          if sum(index_check_up~=0)
             disp('some variables exceed the upper limit')
              output(index_upper_limit(index_check_up)) =in_domain(index_upper_limit(index_check_up));           
          end
     end
        
   if ~isempty(index_lower_limit)
       index_check_low= find(x(index_lower_limit)<lower_limit(index_lower_limit));
       if sum(index_check_low~=0)
           disp('some variables exceed the lower limit')
              output(index_lower_limit(index_check_low)) =in_domain(index_lower_limit(index_check_low));
       end         
   end
   
 else % we backtrack (keep the variables) inside the domain
       if ~isempty(index_upper_limit)
       index_check_up= find(x(index_upper_limit)+Delta_x(index_upper_limit)>upper_limit(index_upper_limit));
         if sum(index_check_up~=0)
           %                s_safe =  .99*min(x(index_check_up)./Delta_x(index_check_up)); % the mininum s that guarantee lambda remain non negative  
           %   s      = (min(1,s_safe));
               Delta_x(index_upper_limit(index_check_up)) = ...
                upper_limit(index_upper_limit(index_check_up))-abs((1-theta)*x(index_upper_limit(index_check_up)))...
                 -x(index_upper_limit(index_check_up));
         end
       end
     
      
      if ~isempty(index_lower_limit)
         index_check_low=find(x(index_lower_limit)+Delta_x(index_lower_limit)<lower_limit(index_lower_limit));
         if sum(index_check_low~=0)
%                 s_safe =  .99*min(-x(index_check_low)./Delta_x(index_check_low)); % the mininum s that guarantee lambda remain non negative  
%                 s      = (min(1,s_safe));
                Delta_x(index_lower_limit(index_check_low)) = ...
                 lower_limit(index_lower_limit(index_check_low))+abs((1-theta)*x(index_lower_limit(index_check_low)))...
                  -x(index_lower_limit(index_check_low))  ;
       end         
      end
      
      output= Delta_x;
end