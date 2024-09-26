function f = sym2func(symbolic_Array,~,input_arguments,option)
% sym2func(symbolic_Array,not_important,input_arguments,option)
% this function convert from symbolic variable to handle function.
% It does the job of matlabFunction but faster and works only when all 
% input arguments are scalar values
% symbolic_Array : the matrix or vector we want to convert it to handle
% function
% not_important: accept any thing it does not used in the function
% options.
% 1- option.function_structure: '@(in)[matrix]'(default), '@(in){cell}',
%'@(in1,in2)[matrix]', '@(in1,in2){cell}'. This how we convert the
% symbolic variables to handle function
% 2- option.blockPartitionSize : used only with '@(in){cell}' or '@(in1,in2){cell}'.
% this option presents gives the block dimensions. For example: 
% symbolic_Array =[x y z ; 1  2 x+y] and option.blockPartitionSize= [ 1 2]
% the partitioned cells are c{1,1}=[ x y]  c{1,2} = z   c{2,1}= [ 1 2]  c{2,2}=x+y 



if isempty(input_arguments)
    error('No input arguments')
end

if ~isvector(input_arguments)
    error('Input arguments must be vector')
end


if  strcmp(class(symbolic_Array),'double')&&~isempty(symbolic_Array) % to deal with double
symbolic_Array=simplify(cell2sym(num2cell(symbolic_Array)));
end


for options_check=1 
 if exist('option')
    if isfield(option,'function_structure'); %
    function_structure = option.function_structure;  
    end
    if isfield(option,'blockPartitionSize'); %
    r_max = option.blockPartitionSize(1); 
        c_max = option.blockPartitionSize(2);  

    end
  end

    if  ~exist('r_max')
       r_max=max(size(symbolic_Array));
       c_max = r_max;  

    end

   if ~exist('function_structure')
      if iscell(input_arguments)
        function_structure='@(in)[matrix]';
      else
       function_structure= '@(in1,in2)[matrix]';
      end
   else 
       switch function_structure
       case {'@(in)[matrix]','@(in){cell}'}
            if  ~iscell(input_arguments)
                input_arguments={input_arguments};
            end
       end
                
       
   end
 end

switch function_structure
    case {'@(in)[matrix]','@(in){cell}'}
         n= length(input_arguments{:});% number of the input arguments;
         n_sym= sum(cellfun(@(cell) strcmp(class(cell),'char'),input_arguments));% chech if all inputs in the cell are symbolic variables
         if n_sym~=0
         error('The input arguments must be symbolic vairables')
         end         
         input_arguments=cell2sym(input_arguments);
       
        
        if ~isempty(symbolic_Array)   
            %%%%%%%%%%% convert the inputs name into input
            input= sym('input',[n 1],'real');  
            for i = 1: n
            symbolic_Array   =subs(symbolic_Array,input_arguments(i),input(i));
            end
           %%%%%%%%%%%
           switch function_structure
            case '@(in)[matrix]'               
              symbolic_Array = sym2string_matrix(symbolic_Array);
              symbolic_Array =inputi2input_i_(symbolic_Array,n);

            case '@(in){cell}'   
                 symbolic_Array=partition(symbolic_Array,r_max,c_max); %%symbolic_Array is cells of symbolic variables blocks
               %%%%% convert all symbolic variables blocks to string blocks
               symbolic_Array= cellfunc(@(cell)sym2string_matrix(symbolic_Array),symbolic_Array);
               symbolic_Array= cellfunc(@(cell)inputi2input_i_(symbolic_Array),symbolic_Array) ;
            end
        
        else 
        symbolic_Array =  'zeros(0,0)';
        f = str2func("@(input)"+ symbolic_Array); 
        return
        end
        %%%% returun part
         switch function_structure
            case '@(in)[matrix]'               
              f = str2func("@(input)"+ symbolic_Array); 
            case '@(in){cell}'   
              f = cellfun(@(cell)"@(input)"+cell,symbolic_Array); % here we have handle functions in each cell; 

         end
       
        
             
    case {'@(in1,in2)[matrix]','@(in1,in2){cell}'}
        
        if ~strcmp(class(input_arguments),'sym')% only we accept symbolic variables
           error('The input arguments must be symbolic vairables')
        end
         %%%%%% make the input argument as string 
         if ~isrow(input_arguments)
            input_arguments =sym2cell(input_arguments)';
            input_arguments=cell2sym(input_arguments);
         end
         
        
          input_arguments= sym2string_matrix(input_arguments);
          input_arguments=char(input_arguments);
        
         if length(input_arguments)==1
            input_arguments= strcat('@(',input_arguments,')'); 
         else
           input_arguments(1)='(';
           input_arguments(end)=')';
           input_arguments = "@"+input_arguments;
         end
         %%%%%%%%%%%%%
         
          if ~isempty(symbolic_Array)   
            switch function_structure
             case '@(in1,in2)[matrix]'               
              symbolic_Array = sym2string_matrix(symbolic_Array);
             case '@(in1,in2){cell}'   
                 symbolic_Array=partition(symbolic_Array,r_max,c_max); %%symbolic_Array is cells of symbolic variables blocks
                 %%%%% convert all symbolic variables blocks to string blocks
                 symbolic_Array= cellfunc(@(cell)sym2string_matrix(symbolic_Array),symbolic_Array) ;
                 
            end
        
          else  
              symbolic_Array =  'zeros(0,0)';
               f = str2func(input_arguments+ symbolic_Array); 
               return
          end
        %%%% returun part
         switch function_structure
            case '@(in1,in2)[matrix]'               
              f = str2func(input_arguments+ symbolic_Array); 
            case '@(in1,in2){cell}' 
              f = cellfun(@(cell)input_arguments+cell,symbolic_Array); % here we have handle functions in each cell; 
         end
             
        
    otherwise
      disp("function_structure:'@(in)[matrix]', '@(in1,in2)[matrix]', '@(in){cell}',or '@(in){cell}' ");
      return
end

   
 
end


function parts     = partition(matrix,r_size,c_size)
% devide the matrix into cells each has specific size
% r_size : is the maximum number of rows in each cell
% c_size : is the maximum number of rows in each c_size
 [row col]=size(matrix); %number of row and column in the matrix
 
 row_index_low =1:r_size:row;
 row_index_up =r_size:r_size:row;
 if length(row_index_low)> length(row_index_up)
 row_index_up(end+1) = row;
 end
 
 col_index_low =1:c_size:col;
 col_index_up =c_size:c_size:col;
 if length(col_index_low)> length(col_index_up)
 col_index_up(end+1) = col;
 end

 
 parts=cell(length(row_index_low),length(col_index_up));
for i=1:length(row_index_low)
      for j=1:length(col_index_up)
          parts{i,j}= matrix(row_index_low(i):row_index_up(i),col_index_low(j):col_index_up(j));               
      end
end


end 


function symbolic_Array   =sym2string_matrix(symbolic_Array);
empty_string     =""; 
symbolic_Array   =char(symbolic_Array);
symbolic_Array   =strrep(symbolic_Array,'matrix([','');
symbolic_Array   =strrep(symbolic_Array,'], [',';');
symbolic_Array   =strrep(symbolic_Array,']])',']');
symbolic_Array   =empty_string+symbolic_Array;
end

function symbolic_Array   =inputi2input_i_(symbolic_Array,n);
% n is number of the inputs
for i=n:-1:1
    symbolic_Array=strrep(symbolic_Array,"input"+i,"input("+i+")");
end
end

