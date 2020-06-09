function algorithm = algorithm_check(algorithm_name)
   
  if   strcmp(algorithm_name,'slack-barrier')
       algorithm_name  ='primal-dual';
  end
 algorithm=algorithm_name;
end
