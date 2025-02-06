%% DEMO_0005_Multi_Morph_Cylindrical_Gyroid_
% This is a demo for:
% 
% * Building geometry for multi-morphology lattices of different gyroid
% structures in cylindrical domain using hybrid formulation.
%
% # Example-1: Utilizes hybrid formulation in axial direction.
% # Example-2: Utilizes hybrid formulation in circumferential direction.
%
%%
%
%  Change log:
%  2023/11/15 MV Created  
%  2024/02/1 MV Sorted for publishing
%  2024/11/21 MV updated with multiMorph
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
res=100; %Resolution

L=9; %Length size
R=3; %Radius size

kappa =8; % kappa controls the lengthscale of transition between lattices

%% Example-1: Axial Transition (Figure-5(a))

inputStruct_A.L=L; % characteristic length
inputStruct_A.R=R;
inputStruct_A.Ns=res; % number of sampling points
inputStruct_A.isocap=1; %Option to cap the isosurface
inputStruct_A.surfaceCase='g'; %Surface type

inputStruct_B = inputStruct_A;
inputStruct_C = inputStruct_A;

% Set parameters for individual gyroid
inputStruct_A.numPeriods=[8 8 8]; %Number of periods in each direction
inputStruct_A.levelset=-0.7; %Isosurface level
inputStruct_A.gradiantF=0; %Gradiant Factor

inputStruct_B.numPeriods=[5 5 5];
inputStruct_B.levelset=-0.8; 
inputStruct_B.gradiantF=0; %Gradiant Factor

inputStruct_C.numPeriods=[6 6 6];
inputStruct_C.levelset= -0.6; 
inputStruct_C.gradiantF=0 ; %Gradiant Factor 

% Compute individual gyroids
% No need to store faces and vertices, only require underlying S,
% grid coordinates, and levelset values
[F,V,C,S_A,X,Y,Z,~,~]=CylindricalTPMS(inputStruct_A);
[~,~,~,S_B,~,~,~,~,~]=CylindricalTPMS(inputStruct_B);
[~,~,~,S_C,~,~,~,~,~]=CylindricalTPMS(inputStruct_C);

% Define the central location of each individual gyroids in space
% E.g., At center_A, the gyroid will definitely correspond to input_A.
% As we move away from center_A, it will slowly transition into other
% gyroids with input_B and input_C. 
center_A = [0, 0, L/6];
center_B = [0, 0, 3*L/6];
center_C = [0, 0, 5*L/6];

% Calculating the multi-morphology surface
Data.S = {S_A; S_B; S_C};
Data.C = {center_A; center_B; center_C};
Data.L = {inputStruct_A.levelset; inputStruct_B.levelset; inputStruct_C.levelset};
Data.K = kappa;

graded_S = GRFBmultiMorph(X,Y,Z,Data);

% Compue isosurface
graded_levelset = 0;
[f,v] = isosurface(X,Y,Z,graded_S,graded_levelset);

% Compute isocaps
[fc,vc] = isocaps(X,Y,Z,graded_S,graded_levelset,'enclose','below');

% Join, merge, and clean unused
[f,v,c] = FV_arrange(f,v,fc,vc);

% Visualize
Hybrid_vizualize(f,v,c,[], Data.C);

%% Example-2: Circumferential Transition (Figure-5(b))

inputStruct_A.L=L; % characteristic length
inputStruct_A.R=R;
inputStruct_A.Ns=res; % number of sampling points
inputStruct_A.isocap=1; %Option to cap the isosurface
inputStruct_A.surfaceCase='g'; %Surface type

inputStruct_B = inputStruct_A;
inputStruct_C = inputStruct_A;

% Set parameters for individual gyroid
inputStruct_A.numPeriods=[8 8 8]; %Number of periods in each direction
inputStruct_A.levelset=-0.7; %Isosurface level
inputStruct_A.gradiantF=0; %Gradiant Factor
levelset_A=inputStruct_A.levelset; 

inputStruct_B.numPeriods=[5 5 5];
inputStruct_B.levelset=-0.8; 
inputStruct_B.gradiantF=0; %Gradiant Factor
levelset_B=inputStruct_B.levelset; 

inputStruct_C.numPeriods=[6 6 6];
inputStruct_C.levelset= -0.6; 
inputStruct_C.gradiantF=0 ; %Gradiant Factor
levelset_C=inputStruct_C.levelset;

% Compute individual spinodoids
% No need to store faces and vertices, only require underlying S,
% grid coordinates, and levelset values
[F,V,C,S_A,X,Y,Z,r,theta]=CylindricalTPMS(inputStruct_A);
[~,~,~,S_B,~,~,~,~,~]=CylindricalTPMS(inputStruct_B);
[~,~,~,S_C,~,~,~,~,~]=CylindricalTPMS(inputStruct_C);

% Define the central location of each individual gyroids in space
% E.g., At center_A, the gyroid will definitely correspond to input_A.
% As we move away from center_A, it will slowly transition into other
% gyroids with input_B and input_C. 

% Converting the grid to cylindrical coordinates
theta_A = 0;
theta_B = 2*pi/3; 
theta_C = 4*pi/3; 

center_A = [R*cos(theta_A), R*sin(theta_A), L/2];
center_B = [R*cos(theta_B), R*sin(theta_B), L/2];
center_C = [R*cos(theta_C), R*sin(theta_C), L/2];

% Calculating the multi-morphology surface
Data.S = {S_A; S_B; S_C};
Data.C = {center_A; center_B; center_C};
Data.L = {inputStruct_A.levelset; inputStruct_B.levelset; inputStruct_C.levelset};
Data.K = kappa;
Data.T = 2; % a circumferential transition

graded_S = GRFBmultiMorph(X,Y,Z,Data);

% Compue isosurface
graded_levelset = 0;
[f,v] = isosurface(X,Y,Z,graded_S,graded_levelset);

% Compute isocaps
[fc,vc] = isocaps(X,Y,Z,graded_S,graded_levelset,'enclose','below');

% Join, merge, and clean unused
[f,v,c] = FV_arrange(f,v,fc,vc);

% Visualize
Hybrid_vizualize(f,v,c,[], Data.C);

%% 
% _*LatticeWorks footer text*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
% 
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
