function Install_IPsolver()
% Run this file to install the IPsolver


if ispc
   slash_type = '\';
else
   slash_type = '/';
end

IPsolverpath = mfilename('fullpath');
slashes = strfind(IPsolverpath, slash_type);
IPsolverpath = IPsolverpath(1:slashes(end));
IPsolverpath = strcat(IPsolverpath,'functions',slash_type) 
addpath(genpath(IPsolverpath))

disp('The directory')
disp(['    ' IPsolverpath])      
disp('has been added to the MATLAB path.')
   savepath
% disp('To use IPsolver consistently, save this new path definition');
% disp('To do this, type the command' );
% disp('savepath' );
end

