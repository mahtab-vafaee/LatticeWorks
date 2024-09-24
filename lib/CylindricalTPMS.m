function [varargout] = CylindricalTPMS(inputStruct)

% function [varargout]=CylindricalTPMS(inputStruct);
%
% -----------------------------------------------------------------------
% This function is designed for generating TPMS lattices in cylindrical
% coordinates. 
% Grid data, faces, vertices, face color data, generated surfaces, and 
% cylindrical coordinates of the grids can be produced as output of this
% function.
%
% _*Name*_
%
% License: <hyperlink to license>
% Author: _Mahtab Vafaee_, <mahtab.vafaee@gmail.com>
%
%  Change log:
%  2024/1/31 MV Created
%  2024/2/1 MV Fixed the variable input cases
% -----------------------------------------------------------------------
%% Parse input

%Create default structure
defaultInputStruct.L=1; % characteristic length
defaultInputStruct.R=1; % characteristic radious
defaultInputStruct.Ns=80; % number of sampling points
defaultInputStruct.isocap=1; %Option to cap the isosurface
defaultInputStruct.thetaMax=2*pi; % characterictic angle
defaultInputStruct.surfaceCase='g';
defaultInputStruct.numPeriods=[1 1 1]; 
defaultInputStruct.levelset=0.5;
defaultInputStruct.surfaceSide=1;
defaultInputStruct.phaseShift=0*[pi pi pi]; 
defaultInputStruct.gradiantF=0;

%Complete input with default if incomplete
[inputStruct]=structComplete(inputStruct,defaultInputStruct,1); %Complement provided with default if missing or empty

%Get parameters from input structure
L = inputStruct.L; % characteristic length
R = inputStruct.R;  % radius
Ns = inputStruct.Ns; % number of sampling points
isocap= inputStruct.isocap; 
surfaceCase=inputStruct.surfaceCase; 
numPeriods=inputStruct.numPeriods; 
levelset=inputStruct.levelset;
phaseShift=inputStruct.phaseShift; 
GF=inputStruct.gradiantF;

%% Generation of points

%Define coordinate limits
rMin=0+phaseShift(1);
rMax=rMin+2*pi;
thetMin=0+phaseShift(2);
thetMax=thetMin+inputStruct.thetaMax;
zMin=0+phaseShift(3);
zMax=zMin+2*pi;

%Create coordinates
if (length(L) == 1)
    xRange=linspace(0,rMax,Ns);
    yRange=linspace(0,thetMax,Ns);
    zRange=linspace(0,zMax,Ns);
    [r,theta,Z]=meshgrid(xRange,yRange,zRange);
end

% Convert cylindrical coordinates to Cartesian coordinates
X = r .* cos(theta);
Y = r .* sin(theta);
Z = Z ;

%Scaling coordinates
if (length(L) == 1)
    X=((X./abs(rMax-rMin)).*R);
    Y=((Y./abs(rMax-rMin)).*R);
    Z=((Z./abs(zMax-zMin)).*L);

else
    X=((X./abs(rMax-rMin)).*L(1));
    Y=((Y./abs(thetMax-thetMin)).*L(2));
    Z=((Z./abs(zMax-zMin)).*L(3));

end

%% Generation of level sets and exporting
S=(sin(X.* numPeriods(1)).*cos(Y.* numPeriods(2)))...
    +(sin(Y.* numPeriods(2)).*cos(Z* numPeriods(3)))...
    +(cos(X.* numPeriods(1)).*sin(Z.* numPeriods(3))); 

% Apply grading
S = S.*(1 + (-GF*(Z/zMax))); %Linear gradient in z-direction

%Iso-surface
switch inputStruct.surfaceSide
    case 0
        [f,v,c]=getSurface(X,Y,Z,S,levelset);
        if isocap==1
            [fc1,vc1,cc1]=getCaps(X,Y,Z,S,levelset);
            [fc2,vc2,cc2]=getCaps(X,Y,Z,-S,-levelset);
            [f,v,c]=joinElementSets({f,fc1,fc2},{v,vc1,vc2},{c,cc1,cc2+6});
        end
    case 1 
        [f,v,c]=getSurface(X,Y,Z,S,levelset);        
        if isocap==1
            [fc,vc,cc]=getCaps(X,Y,Z,S,levelset);
            [f,v,c]=joinElementSets({f,fc},{v,vc},{c,cc}); %Join sets
        end
    case -1
        [f,v,c]=getSurface(X,Y,Z,-S,-levelset);        
        if isocap==1 
            [fc,vc,cc]=getCaps(X,Y,Z,-S,-levelset);
            [f,v,c]=joinElementSets({f,fc},{v,vc},{c,cc}); %Join sets
        end        
end

%% Merge nodes and clean-up mesh 

%Merge nodes
[f,v]=mergeVertices(f,v); 

%Check for unique faces
[~,indUni,~]=unique(sort(f,2),'rows');
f=f(indUni,:); %Keep unique faces
c=c(indUni);

%Remove collapsed faces
[f,logicKeep]=patchRemoveCollapsed(f); 
c=c(logicKeep);

%Remove unused points
[f,v]=patchCleanUnused(f,v); 

%Invert faces
f=fliplr(f); 

%% Collect output
varargout{1}=f;
varargout{2}=v;
varargout{3}=c;
varargout{4}=S;
varargout{5}=X;
varargout{6}=Y;
varargout{7}=Z;
varargout{8}=r;
varargout{9}=theta;

end
%% 
% _*LatticeWorks footer text*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
% 
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
