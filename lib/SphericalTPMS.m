function [varargout]=SphericalTPMS(inputStruct)

% function [varargout]=SphericalTPMS(SphericalTPMS);
%
% -----------------------------------------------------------------------
% This function is used to generate TPMS structures in spherical
% coordinates.
%
% _*Name*_ 
% 
% License: <hyperlink to license>
% Author: _Mahtab Vafaee_, <mahtab.vafaee@gmail.com>
%
%  Change log:
%  2024/1/31 MV Created  
% -----------------------------------------------------------------------
%% Parse input
%Create default structure
defaultInputStruct.L=1; % characteristic length
defaultInputStruct.Ns=80; % number of sampling points
defaultInputStruct.isocap=1; %Option to cap the isosurface
defaultInputStruct.surfaceCase='g';
defaultInputStruct.numPeriods=[1 1 1]; 
defaultInputStruct.levelset=0.5;
defaultInputStruct.surfaceSide=1;
defaultInputStruct.phaseShift=0.*[pi pi pi]; 
defaultInputStruct.GF=0; 

%Complete input with default if incomplete
[inputStruct]=structComplete(inputStruct,defaultInputStruct,1); %Complement provided with default if missing or empty

%Get parameters from input structure
L = inputStruct.L; % characteristic length
Ns = inputStruct.Ns; % number of sampling points
isocap= inputStruct.isocap; 
surfaceCase=inputStruct.surfaceCase; 
numPeriods=inputStruct.numPeriods; 
levelset=inputStruct.levelset;
phaseShift=inputStruct.phaseShift; 

if numel(numPeriods)==1
    numPeriods=numPeriods*ones(1,3);
end

%% Generation of points
%Define coordinate limits
xMin=0+phaseShift(1);
xMax=xMin+2*pi*numPeriods(1);
yMin=0+phaseShift(2);
yMax=yMin+2*pi*numPeriods(2);
zMin=0+phaseShift(3);
zMax=zMin+2*pi*numPeriods(3);

%Create coordinates
xRange=linspace(xMin,xMax,Ns);
yRange=linspace(yMin,yMax,Ns);
zRange=linspace(zMin,zMax,Ns);
[X,Y,Z]=meshgrid(xRange,yRange,zRange);

%Calculate 3D image data
S=triplyPeriodicMinimal(X(:),Y(:),Z(:),surfaceCase);        
S=reshape(S,size(X));

%Scaling coordinates
switch length(L)
    case 1
        X=(((X./abs(xMax-xMin))-1/2).*L);
        Y=(((Y./abs(yMax-yMin))-1/2).*L);
        Z=(((Z./abs(zMax-zMin))-1/2).*L);

    case 3
        X=(((X./abs(xMax-xMin))-1/2).*L(1,1));
        Y=(((Y./abs(yMax-yMin))-1/2).*L(1,2));
        Z=(((Z./abs(zMax-zMin))-1/2).*L(1,3));
end

%% Gradient surface
switch length(L)
    case 1
        S=S.*(1-inputStruct.GF*Z/L);
    case 3
        S=S.*(1-inputStruct.GF*Z/L(1,3));
end


%% Generation of level sets and exporting
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

end
%% 
% _*LatticeWorks footer text*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
% 
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
