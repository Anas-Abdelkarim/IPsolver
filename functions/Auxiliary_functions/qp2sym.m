function [solution] = qp2sym(H,f,A,b,Aeq,beq,lb,ub,X0,option); 
%qp2sym(H,f,A,b,Aeq,beq,LB,UB,X0,options); 
%min .5*x'Hx+f'x  s.t  Ax<=b; Aeq<=beq ;LB<=x<UB

if isempty(H) 
   [n, ~] =  size(f);  %n: number of variales
    H=zeros(n,n);
else 
     [n, ~] =  size(H); 
end
if isempty(f)
    f= zeros(n,1);
else
   if isrow(f)
        f=f';
   end
end
if isempty(Aeq)
 Aeq=double.empty(0 ,n);
 beq=double.empty(0 ,1); 
end

if isempty(A)
 A=double.empty(0 ,n);
 b=double.empty(0 ,1); 
end
if ~isempty(lb)
  f_i_1=[lb<=x];
else 
  f_i_1=[];
end

if ~isempty(lb)
  f_i_2=[x<=ub];
else 
  f_i_2=[];
end


x = sym('x' ,[n,1],'real');
decision_variables= x;
f_0= .5*(x'*H*x)+f'*x;
equality= Aeq*x==beq;
f_i= [ A*x<=b; f_i_1;f_i_2];
parameters= [];
parameters_subs=[];
warm_start =[]; 


algorithm= option.algorithm;
KKT_matrix = KKT(decision_variables,f_0,f_i,equality,parameters,algorithm,option);
solution   = IPsolver(KKT_matrix,X0,parameters_subs,warm_start);
end
