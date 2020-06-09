function Uninstall_IPsolver()
% Run this file to uninstall the IPsolver


if ispc
   slash_type = '\';
else
   slash_type = '/';
end

IPsolverpath = mfilename('fullpath');
slashes = strfind(IPsolverpath, slash_type);
IPsolverpath = IPsolverpath(1:slashes(end));
IPsolverpath = strcat(IPsolverpath,'functions',slash_type) 
rmpath(genpath(IPsolverpath))

disp('The directory')
disp(['    ' IPsolverpath])
 disp('has been removed form the MATLAB path.')
savepath
end

