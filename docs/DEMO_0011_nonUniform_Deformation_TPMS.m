%% DEMO_0011_nonUniform_Deformation_TPMS
% This is a demo for:
% 
% * Building geometry for an arbitrary non-uniform deformation matrix
% applied on gyroid structure.
%
% This demo contains:
%
% # Case-1: Gyroid lattice under twist deformation.
% # Case-2: Gyroid lattice under rotation deformation.
% # Case-3: Gyroid lattice under bending deformation.
% # Case-4: Gyroid lattice under shear deformation.
% # Case-5: Gyroid lattice under stretch deformation.
%
%%
%
%  Change log:
%  2023/11/15 MV Created  
%  2024/02/2 MV  Edited
%  2024/11/19 MV Edited
%  2025/02/06 MV Added lib path
% ----------------------------------------------------------------------
%%

clear; close all; clc;

%% Adding lib path so functions are known 

mainPath=fileparts(mfilename('fullpath')); %Get the  path
addpath(fullfile(fileparts(mainPath),'lib')); %Add lib path 
addpath(fullfile(fileparts(mainPath),'lib_ext')); %Add external lib path 

%%
% Plot settings
cMap=jet(250);
faceAlpha1=1;
faceAlpha2=0.65;
edgeColor1='none';
edgeColor2='none';
fontSize=20; 
pColors=gjet(6);

%% Control parameters

n= 100; % resolution 
boxDim = [1, 1, 2]; % dimenasion 
l=0.9;% levleset

% Define deformation control parameters
% Rotation & Twist
bendAngle = pi/2; 
rInner = 2;
barHeight = bendAngle*rInner; 

% Radius of curvature for bending
bend_radius = 40; 

% Shear factor (controls the amount of shear)
shear_factor = 0.2; 

% Stretch factor (controls the amount of stretch)
stretch_factor = 0.2; 
necking_factor = 0.3; % Controls the necking effect

%% Create an origiional grid

DefType='Bending'; % {'Twist','Rotate', 'Bending', 'Shear', 'Stretch'}

switch DefType
    case 'Rotate'
        xRange = linspace(-pi,0,boxDim(1,1)*n)-rInner;
        yRange = linspace(0,2*pi,boxDim(1,2)*n);
        zRange = linspace(0,barHeight,boxDim(1,3)*n);

    case {'Twist', 'Bending', 'Shear'}
        xRange = linspace(-pi,pi,boxDim(1,1)*n);
        yRange = linspace(-pi,pi,boxDim(1,1)*n);
        zRange = linspace(0,10,boxDim(1,1)*n);

    case 'Stretch'
        xRange = linspace(-pi,pi,boxDim(1,1)*n);
        yRange = linspace(-pi,pi,boxDim(1,1)*n);
        zRange = linspace(-5,5,boxDim(1,1)*n);
end

[X,Y,Z]=meshgrid(xRange,yRange,zRange);

%% Applying deformation

%Deformation matrix
a=linspace(0,bendAngle,size(Z,3));
Xp=X; Yp=Y; Zp=Z;

switch DefType
    case 'Rotate'
        for q=1:1:size(Z,3)
            R = euler2DCM([0 a(q) 0]);
            x = X(:,:,q);
            y = Y(:,:,q);
            z = zeros(size(x));
            v = [x(:) y(:) z(:)];
            %     vp = v*R;
            vp = (R*v')';

            xp = reshape(vp(:,1),size(x));
            yp = reshape(vp(:,2),size(x));
            zp = reshape(vp(:,3),size(x));

            Xp(:,:,q)=xp;
            Yp(:,:,q)=yp;
            Zp(:,:,q)=zp;
        end

    case 'Twist'
        for q=1:1:size(Z,3)
            R = euler2DCM([0 0 a(q)]);
            x = X(:,:,q);
            y = Y(:,:,q);
            z = Z(:,:,q);
            v = [x(:) y(:) z(:)];
            %     vp = v*R;
            vp = (R*v')';

            xp = reshape(vp(:,1),size(x));
            yp = reshape(vp(:,2),size(x));
            zp = reshape(vp(:,3),size(x));

            Xp(:,:,q)=xp;
            Yp(:,:,q)=yp;
            Zp(:,:,q)=zp;
        end

    case 'Bending'
        % Angular deformation based on Z coordinates
        theta = Z./bend_radius; 

        % Apply rotation the grid in xz surface around the Y-axis
        Xp = X.*cos(theta)+Z.*sin(theta);
        Zp = -X.*sin(theta)+Z.*cos(theta);
        Yp = Y; % Y remains unchanged
        
    case 'Shear'
        % Shear in the X-direction, normal to Z
        Xp = X - shear_factor*Z; 
        Yp = Y; 
        Zp = Z; 

    case 'Stretch'
        % Stretching in Z-direction;
        Zp = Z+stretch_factor*Z; 

        neck_profile = 1 - necking_factor * exp(-((Z / max(abs(zRange))).^2)); % Gaussian necking
        % Applying necking profile to X and Y
        Xp = X .* neck_profile; 
        Yp = Y .* neck_profile; 
        
end

%% Create deformed grid

Vo = [X(:) Y(:) Z(:)]; %Origional grids
Vp = [Xp(:) Yp(:) Zp(:)]; %Deformed grids

%% Evaluate triply periodic function

% calculate the grid value to create gyroid field
freq1=5; %period number
S=(sin(freq1.*X).*cos(freq1.*Y))+...
    (sin(freq1.*Y).*cos(freq1.*Z))+...
    (cos(freq1.*X).*sin(freq1.*Z)); 

% Construct deformed iso-surface
Sn=S;
[Fi,Vi] = isosurface(Xp,Yp,Zp,Sn,l);
[Fc,Vc] = isocaps(Xp,Yp,Zp,Sn,l);
[F,V] = joinElementSets({Fi,Fc},{Vi,Vc});
[F,V] = mergeVertices(F,V);
F=fliplr(F);

% Construct origional iso-surface
[Fi,Vi] = isosurface(X,Y,Z,Sn,l);
[Fc,Vc] = isocaps(X,Y,Z,Sn,l);
[Fj,Vj] = joinElementSets({Fi,Fc},{Vi,Vc});
[Fj,Vj] = mergeVertices(Fj,Vj);
Fj=fliplr(Fj);

%% Visualize deformed surface 
cFigure;

subplot(1,3,1);hold on;
title('Original Structure','FontSize',fontSize);
gpatch(Fj,Vj,'kw','none',1);
axisGeom(gca, fontSize);
camlight headlight;

hp1=subplot(1,3,2);hold on;
title('Deformed Grids','FontSize',fontSize);
plotV(Vo,'k.','MarkerSize',30);
plotV(Vp,'r.','MarkerSize',30);
axisGeom(gca, fontSize);
camlight headlight;

hp2=subplot(1,3,3);hold on;
title('Deformed Gyroid','FontSize',fontSize);
gpatch(F,V,'kw','none',1);
axisGeom(gca, fontSize);
camlight headlight;
drawnow; 

%% 
% _*LatticeWorks footer text*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
% 
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
