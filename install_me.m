
%% Add library paths so functions are known

mainPath=fileparts(mfilename('fullpath')); %Get the main path
addpath(fullfile(mainPath,'lib')); %Add lib folder containing custom functions
addpath(fullfile(mainPath,'lib_ext')); %Add lib_ext folder containing external/3rd party functions
savepath; % Saving the newly added paths (might only work if user has admin rights)
