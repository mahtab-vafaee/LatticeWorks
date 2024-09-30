clear; close all; clc;

%%

toolboxPath=fileparts(fileparts(fileparts(mfilename('fullpath'))));
docPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'docs');
libPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'lib');
docToolsPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'docs','docTools');

pathNames={docToolsPath,libPath,docPath};

%% Prepare boilerplate

licenseBoilerPlate=fullfile(toolboxPath,'codeFooter.txt');
[T]=txtfile2cell(licenseBoilerPlate);

footerTargetText='% _*LatticeWorks footer text*_ ';

%%

licenseLink='https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE';

%Add comment symbols
for q=1:1:numel(T)
    T{q}=['% ',T{q}];
end

%Add target header
T=[{'%% '};{footerTargetText};{'% '};{['% License: <',licenseLink,'>']};T(1:end)];

%%

fileExtension='.m';

%%

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
        T_now(end+1:end+numel(T))=T;
        cell2txtfile(fileName,T_now,0);        
    end
end

%% 
% _*LatticeWorks footer text*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
% 
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
