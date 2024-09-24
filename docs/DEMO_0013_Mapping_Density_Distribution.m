%% DEMO_0013_Mapping_Density_Distribution
% This is a demo for:
% 
% * Building a graded lattice structure, to map a specific 
% structural properties on a rectangular domain, 
% e.g. in this demo, it is mapping a density distribution field.
%
%%
%
%  Change log:
%  2023/11/15 MV Created  
%  2024/02/06 MV Edited
%  2024/04/15 MV/KM Revised
% ----------------------------------------------------------------------
%%
clc; clear; close all;

%% Plot settings

cMap=jet(250);
faceAlpha1=0.5;
faceAlpha2=0.65;
edgeColor1='none';
edgeColor2='none';
fontSize=25; 
pColors=gjet(6);

%% Using topology optimisation to create a density field
% define box size, volume fraction as inputs for top.m
s=2; 
nelx = s*64;
nely = s*24;
nelz = s*24;

volfrac = 0.5;
penal = 3;
rmin = 5;

%% TOP: calculating the density distribution map
[rho0] = top(nelx,nely,volfrac,penal,rmin);

%% extending through the z-direction
ext = [1, 1, nelz];
rho = repmat(rho0, ext); 

%% Visualize 3D model and density data
sv3(rho); colormap warmcold; 

%% Defining desired spacing (e.g. as per desired spatial frequency) for grids

spatialWaveLength = 5; % Nominal spatial wave length
samplingRate = 3; % Rate per "wave"
v = spatialWaveLength/samplingRate;

%% Creating a grid with desired voxel size for Gyroid image

w=linspace(0,nelx,nelx);
d=linspace(0,nely,nely);
h=linspace(0,nelz,nelz);

[X,Y,Z]=meshgrid(w,d,h); 

% grids with desired resolution v
[XG, YG, ZG] = meshgrid(linspace(0,nelx,v*nelx),...
    linspace(0,nely,v*nely),...
    linspace(0,nely,v*nely));

%% Interpolating data on the grid

rho_VG = griddata(X,Y,Z,rho,XG(:),YG(:),ZG(:),'linear'); %Interpolate on grid
rho_VG = reshape(rho_VG,size(XG));

%% Visualize interpolated density on the grid
sv3(rho_VG); colormap warmcold;
caxis([0,1]);

%% Mapping the density field to the equivalent gyroid levelSet field

l=(rho_VG-0.5)/-(1/3); % Equivalent levelset image for rho
k = (1*pi)/spatialWaveLength; % number of periods

%% Visualizig levelSet data
cFigure;  
scatter3(XG(:),YG(:),ZG(:),20,l(:),'filled');
axis tight; axis equal; 
colorbar;

%% Evaluate triply periodic function

S=(sin(k.*XG).*cos(k.*YG))+(sin(k.*YG).*cos(k.*ZG))+(cos(k.*XG).*sin(k.*ZG)); 

sv3(S); colormap warmcold;

%% Mapping the density field to the equivalent gyroid levelSet field
% Gyroid "spans" -1.5-1.5, shift so minimum is at 1
i=2.5; % Leads to a minimum of 1 later
Sg=S+i; % S=-S-i;
l=l+i; % l=-l-i; % i:i+s
Sg=Sg./l; 

Sg=Sg-1; % bringing back to zero levelset
% Sg=-Sg; % Negative

%visualize gyroid function field
sv3(Sg); colormap warmcold;

%% Construct iso-surface

[Fg,Vg] = isosurface(XG,YG,ZG,Sg,0.5);
[Fgc,Vgc] = isocaps(XG,YG,ZG,Sg,0);

%% Joining the two surfaces

[Fsn,Vsn,Csn]=joinElementSets({Fg,Fgc},{Vg,Vgc});
[Fsn,Vsn]=patchCleanUnused(Fsn,Vsn); %Remove unused nodes
[Fsn,Vsn]=mergeVertices(Fsn,Vsn); %Merge nodes

%% Visualize surface
cFigure; 
gpatch(Fsn,Vsn,'kw','none',1);
axisGeom; colormap spectral; icolorbar; 
camlight headlight; axis on;
drawnow; 

%% Scallin the Coordinates

scalF=max(squeeze(Vsn(:, 1, :)))/60;
Vsn1=Vsn./scalF;

%% STL Export

TR = triangulation(Fsn,Vsn1(:,1),Vsn1(:,2),Vsn1(:,3));
stlwrite(TR,'3D_Printing_Prototye.stl');
%% 
% _*LatticeWorks footer text*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
% 
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
