%% DEMO_0009_Multi_Morph_Cylindrical_Arrangment
% This is a demo for:
% 
% * Building geometry for multi-morphology TPMS structures (gyroid and diamond)
% in cylindrical arrangement, with transition in different directions.
%
% This demo contains:
%
% # Case-1: TPMS in cylindrical arrangment, radial transition.
% # Case-2: TPMS in cylindrical arrangment, circumferential transition.
% # Case-3: TPMS in cylindrical arrangment, axial transition.
%
%%
%
%  Change log:
%  2023/11/15 MV Created  
%  2024/02/2 MV 
%  2024/11/21 MV Updated with multiMorph
%  2025/02/06 MV Added lib path
% ----------------------------------------------------------------------
%%

clear; close all; clc;

%% Adding lib path so functions are known 

mainPath=fileparts(mfilename('fullpath')); %Get the  path
addpath(fullfile(fileparts(mainPath),'lib')); %Add lib path 
addpath(fullfile(fileparts(mainPath),'lib_ext')); %Add external lib path 

%% Plot settings
fontSize=20;
faceAlpha1=0.8;
markerSize=10;
lineWidth1=3;
lineWidth2=4;
markerSize1=25;

%% Control parameters

res=80; %Resolution

%% Setting-up input parameters for individual lattices

inputStruct_A.size=[0, 3, 2*pi, 3]; % [r_in, r_out, theta, length], for section view :[0, 3, pi, 3]
inputStruct_A.Ns=res; % number of sampling points
inputStruct_A.isocap=1; %Option to cap the isosurface
inputStruct_A.surfaceCase='g'; %Surface type

r1=inputStruct_A.size(1,1); % inner radius 
r2=inputStruct_A.size(1,2); % outter radius
L=inputStruct_A.size(1,4); % height

inputStruct_B = inputStruct_A;
inputStruct_C = inputStruct_A;

% Set parameters for individual gyroid
inputStruct_A.numPeriods=[3 10 2]; %Number of periods in each direction
inputStruct_A.levelset=-0.8 ; %Isosurface level

inputStruct_B.numPeriods=[2 8 2];
inputStruct_B.levelset=-0.9; 
inputStruct_B.surfaceCase='d'; 

inputStruct_C.numPeriods=[2 20 2];
inputStruct_C.levelset=-0.5; 
inputStruct_C.surfaceCase='g'; 

%% Compute individual TPMS
% No need to store faces and vertices, only require underlying S,
% grid coordinates, and levelset values

[F,V,C,S_A,X,Y,Z,r,theta]=TPMS_LCS(inputStruct_A);
[~,~,~,S_B,~,~,~,~,~]=TPMS_LCS(inputStruct_B);
[~,~,~,S_C,~,~,~,~,~]=TPMS_LCS(inputStruct_C);

%% Transition lengthScale and shape
% kappa controls the lengthscale of transition between lattices
kappa =15;

% Transition shape type
transCase = 'c';
switch transCase
    case 'r' % radial
        center_A = [1, 1, 1];
        center_B = [2, 2, 1];
        center_C = [4, 4, 1];

        Data.T = 3;

    case 'c' % circumferential
        theta_A=0;
        theta_B=2/3*pi;
        theta_C=4/3*pi;

        center_A = [r2*cos(theta_A), r2*sin(theta_A), L];
        center_B = [r2*cos(theta_B), r2*sin(theta_B), L];
        center_C = [r2*cos(theta_C), r2*sin(theta_C), L];

    case 'a' % axial
        center_A = [0, 0, 0.5];
        center_B = [0, 0, 1.5];
        center_C = [0, 0, 2.5];
end

%% Calculating the multi-morphology surface

Data.S = {S_A; S_B; S_C};
Data.C = {center_A; center_B; center_C};
Data.L = {inputStruct_A.levelset; inputStruct_B.levelset; inputStruct_C.levelset};
Data.K = kappa;

graded_S = GRFBmultiMorph(X,Y,Z,Data);

%% Compue isosurface
graded_levelset = 0;
[f,v] = isosurface(X,Y,Z,graded_S,graded_levelset);

% Compute isocaps
[fc,vc] = isocaps(X,Y,Z,graded_S,graded_levelset,'enclose','below');

% Join, merge, and clean unused
[f,v,c] = FV_arrange(f,v,fc,vc);

%% Visualize
map=[0.75 0.75 0
    0	0.5 0
    0   0.5 0
    0	0.5	0
    0	0.5	0
    0.75 0.3 0.2
    0.4 0.4 0.9];

Hybrid_vizualize(f,v,c,map,[]);
%% 
% _*LatticeWorks footer text*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
% 
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
