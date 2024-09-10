function f = callfunc(func,input,function_structure)
% f = callfunc(func,input,function_structure))
% this function prepared the output of any function depending on the
% func: any handle function
% input: vector input
% function_structure : '@(in)[matrix]','@(in){cell}','@(in1,in2)[matrix]',
% or  '@(in1,in2){cell}'

% detect the function structure if it is not given
if isempty(func)
    f = [];
    return
end 
if ~exist('function_structure')
function_structure=[];    
end
if isempty(function_structure) 
    if iscell(func)
        if nargin(func{1})==1
            function_structure='@(in){cell}';
        else
            function_structure='@(in1,in2){cell}';
        end               
    else
         if nargin(func)==1
            function_structure='@(in)[matrix]';
        else
            function_structure='@(in1,in2)[matrix]';
        end       
        
        
    end

end



switch function_structure
    case '@(in)[matrix]'
        f = func(input);
    case '@(in1,in2)[matrix]'
        input= num2cell(input);
        f = func(input{:});
    case '@(in){cell}'
        f= [cellfun(@(cell)cell(input),func)]; % here we substitute in each cell
     case '@(in1,in2){cell}'
        input= num2cell(input);
        f= [cellfun(@(cell)cell(input{:}),func)]; % here we substitute in each cell
    otherwise
            disp("function_structure:'@(in)[matrix]', '@(in1,in2)[matrix]', '@(in){cell}',or '@(in){cell}' ")
        f = func(input);
 end
        

    
end
