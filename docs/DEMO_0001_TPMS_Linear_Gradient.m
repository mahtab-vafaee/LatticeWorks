%% DEMO_0001_TPMS_Linear_Gradient
% This is a demo for:
% 
% * Building geometry for linear volume fraction or cell size gradiet TPMS
% lattices

%%
% 
%  Change log:
%  2023/11/15 MV Created  
%  2024/01/29 MV Sorted for publishing
%  2024/09/03 KMM cleaning up, commenting
%  2024/11/21 MV creating more core functions
%  2025/02/06 MV Added lib path
% -----------------------------------------------------------------------

%%

clear; close all; clc;

%% Adding lib path so functions are known 

mainPath=fileparts(mfilename('fullpath')); %Get the  path
addpath(fullfile(fileparts(mainPath),'lib')); %Add lib path 
addpath(fullfile(fileparts(mainPath),'lib_ext')); %Add external lib path 

%% Plot settings

cMap=[0.6*ones(256,2), linspace(1, 0, 256)'];
fontSize=15; 

%% Control parameters
% Creating control parameters depending on the type of demo the user wants
% to run 

gradType= 'levelSet'; % choose between cellSize or levelSet
levelset= 1; %Isosurface level, corresponding to volume fraction

switch gradType
    case 'levelSet' % PAPER, figure 1
        inputStruct.L=[3 1 1]; % characteristic length
        inputStruct.Ns=100; % number of sampling points, resolution
        inputStruct.surfaceCase='g'; %Surface type
        inputStruct.numPeriods=[9 3 3]; %Number of periods in each direction [6 2 2]
        inputStruct.gradType=gradType;
        inputStruct.GF=[-1.2, -0.3]; % Gradient factors

    case 'cellSize' % Coarse version 
        inputStruct.L=[3 1 1]; % characteristic length
        inputStruct.Ns=100; % number of sampling points, resolution
        inputStruct.surfaceCase='g'; %Surface type
        inputStruct.numPeriods=[6 2 2]; %Number of periods in each direction [6 2 2]
        inputStruct.gradType=gradType;
        inputStruct.GF=3; % Gradient factor
        inputStruct.phaseShift=-1/4*pi;
end

%% Create triply periodic minimal surface
[S,X,Y,Z] = gradTPMS(inputStruct);

% Isosurface
[F,V] = isosurface(X,Y,Z,S,levelset);

%Capping ends
[fc,vc]=isocaps(X,Y,Z,S,levelset, 'above'); 

% Join, merge, and clean unused
[f,v,c] = FV_arrange(F,V,fc,vc);

%% Visualize surface

cFigure; 
subplot(1,2,1); hold on;
title('Face labelling');
hp1 = gpatch(f,v,c,'none', 1); 
colormap(gca,gjet(6)); 
axisGeom(gca,fontSize); camlight headlight; 

switch gradType
    case 'levelSet'
        subplot(1,2,2); hold on;        
        title('Gradient');
        hp2 = gpatch(f,v,v(:,1),'none', 1); % color gradient allong X-direction
        hp2.FaceColor = 'interp'; % Use interpolated shading
        colormap(gca,cMap) ; %Set Colormap
        cbh = colorbar;  %Create Colorbar
        cbh.Ticks = linspace(0, 3, 4); %Create 4 ticks 10%-40%
        cbh.TickLabelInterpreter = 'tex';
        cbh.TickLabels = {'10%','20%','30%','40%'} ; %Replace the labels
        % of these 8 ticks with the numbers 1 to 8
    case 'cellSize'        
        subplot(1,2,2); hold on;                
        title('Gradient');
        hp2 = gpatch(f,v,v(:,1),'none', 1);
        hp2.FaceColor = 'interp'; % Use interpolated shading
        colormap(gca,cMap) ; %Set Colormap
        % cbh = colorbar;  %Create Colorbar
end

axisGeom(gca,fontSize); camlight headlight; 
gdrawnow;

%% 
% _*LatticeWorks footer text*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
% 
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
