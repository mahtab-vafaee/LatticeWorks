clear; close all; clc;

%% Path names

toolboxPath=fileparts(fileparts(fileparts(mfilename('fullpath'))));
docPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'docs');
libPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'lib');
docToolsPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'docs','docTools');

pathNames={docToolsPath,libPath,docPath};

%%

fileExtension='.m'; %Extension to remove boilerplate from

footerTargetStart='% _*LatticeWorks footer text*_ '; %Target for removal

%N.B. Removal is up to end from the start so make sure the target is
%appropriately set!!!!!

%% Removing footer text

numPaths=numel(pathNames); 
for q_path=1:1:numPaths
    pathName=pathNames{q_path};     
    files = dir(fullfile(pathName,'*.m'));
    files={files(1:end).name};
    files=sort(files(:));
    numFiles=numel(files);    
    for q_file=1:1:numFiles
        fileName=fullfile(pathName,files{q_file});
        [T_now]=txtfile2cell(fileName);
        targetStartIndex = find(strcmp(footerTargetStart,T_now));       
        if ~isempty(targetStartIndex)            
            targetStartIndex=targetStartIndex(end); %Keep last occurance
            T_now=T_now(1:targetStartIndex-1-1); %NB -1 is used to remove %% above target 
            cell2txtfile(fileName,T_now,0,0);
        end
    end
end

%% 
% _*LatticeWorks footer text*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
% 
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
