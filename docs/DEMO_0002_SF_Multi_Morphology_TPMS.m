%% DEMO_0002_SF_Multi_Morphology_TPMS
% This is a demo for:
% 
% * Building geometry for generating multi-morphology lattice of 
% Gyroid-Diamond in cubic domain using Sigmoid Function (SF)

%%
%
%  Change log:
%  2023/11/15 MV Created  
%  2024/01/31 MV Sorted for publishing
%  2024/11/21 MV Updated with SFmultiMorph
%  2025/02/06 MV Added lib path
% -----------------------------------------------------------------------

%%
clear; close all; clc;

%% Adding lib path so functions are known 

mainPath=fileparts(mfilename('fullpath')); %Get the  path
addpath(fullfile(fileparts(mainPath),'lib')); %Add lib path 
addpath(fullfile(fileparts(mainPath),'lib_ext')); %Add external lib path 

%% Plot settings

fontSize=20;

%% Control parameters

sampleSize=2; %Heigh of the sample
pointSpacing=sampleSize/100; % Resolution

overSampleRatio=1;
numStepsLevelset=ceil(overSampleRatio.*(sampleSize./pointSpacing)); %Number of voxel steps across period for image data (roughly number of points on mesh period)

%% Set-up the input parameters for individual lattices

inputStruct_A.L=[4 2 2]; % characteristic length
inputStruct_A.Ns=numStepsLevelset; % number of sampling points
inputStruct_A.isocap=1; %Option to cap the isosurface
inputStruct_A.surfaceCase='g'; %Surface type

inputStruct_B = inputStruct_A;

% Set parameters for individual lattices 
% Structure-A
inputStruct_A.numPeriods=[5 2 2]; %Number of periods in each direction
inputStruct_A.levelset=-0.1 ; %Isosurface level
inputStruct_A.gradiantF=0 ; %Gradiant Factor within individual structures
levelset_A=inputStruct_A.levelset; 

% Structure-B
inputStruct_B.numPeriods=[6 3 3];
inputStruct_B.levelset=-0.4; 
inputStruct_B.surfaceCase='d'; 
inputStruct_B.gradiantF=0 ; %Gradiant Factor within individual structures
levelset_B=inputStruct_B.levelset; 

%% Compute individual gyroids

% No need to store faces and vertices, only require underlying S,
% grid coordinates, and levelset values
[~,~,~,S_A,X,Y,Z]=triplyPeriodicMinimalSurface(inputStruct_A);
[~,~,~,S_B,~,~,~]=triplyPeriodicMinimalSurface(inputStruct_B);

%% Define the central location of each individual lattices in space
% E.g., At center_A, the structure will definitely correspond to input_A.
% As we move away from center_A, it will slowly transition into other
% structures with input_B. 

center_A = [0.75, 0.5, 0.5];
center_B = [1.25, 0.5, 0.5];

%% Transition lengthscal and shape
% kappa controls the lengthscale of transition between lattices
% Higher kappa => faster transition
% Lower kappa => slower transition
kappa = 5;

%Transition path (shape)
G = X;
G = G/max(G(:));
G = G-(max(G(:))/2);

%% Compute the weitgh functions
Data.S = {S_A; S_B};
Data.C = {center_A; center_B};
Data.L = {inputStruct_A.levelset; inputStruct_B.levelset};
Data.K = kappa;
Data.G = G;

graded_S = SFmultiMorph(Data);

%% Compue isosurface
graded_levelset = 0;
[f,v] = isosurface(X,Y,Z,graded_S,graded_levelset);

% Compute isocaps
[fc,vc] = isocaps(X,Y,Z,graded_S,graded_levelset,'enclose','below');

% Join, merge, and clean unused
[f,v,c] = FV_arrange(f,v,fc,vc);

%% Visualize
Hybrid_vizualize(f,v,c);

%% 
% _*LatticeWorks footer text*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
% 
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
