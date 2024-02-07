function [varargout]= TPMS_LCS(inputStruct)

% function [f,v,c,S,X,Y,Z,r,phi]=TPMS_LCS(inputStruct)
%
% -----------------------------------------------------------------------
% This function is used to generate TPMS lattice in cylindrical
% arranement.
%
% _*Name*_ 
% 
% License: <hyperlink to license>
% Author: _Mahtab Vafaee_, <mahtab.vafaee@gmail.com>
%
%  Change log:
%  2024/1/31 MV Created  
% -----------------------------------------------------------------------
%%
% Parse Input Structure
%Create default structure
defaultInputStruct.r1=1; % characteristic length
defaultInputStruct.r2=2;
defaultInputStruct.theta=pi;
defaultInputStruct.z=2;

defaultInputStruct.a=1; % peridocity in r, theta, z directions
defaultInputStruct.b=1;
defaultInputStruct.c=1;

defaultInputStruct.Ns=80; % number of sampling points
defaultInputStruct.isocap=1; %Option to cap the isosurface
defaultInputStruct.surfaceCase='g';
defaultInputStruct.numPeriods=[1 1 1]; 
defaultInputStruct.levelset=0.5;
defaultInputStruct.surfaceSide=1;
defaultInputStruct.phaseShift=0*[pi pi pi]; 
defaultInputStruct.gradiantF=-4;

%% Complete input with default if incomplete
[inputStruct]=structComplete(inputStruct,defaultInputStruct,1); %Complement provided with default if missing or empty

%Get parameters from input structure
r1= inputStruct.size(1); % characteristic length
r2= inputStruct.size(2);
theta= inputStruct.size(3);
z= inputStruct.size(4);

a=inputStruct.numPeriods(1); % Number of periods in radial direction
b=inputStruct.numPeriods(2); % in theta direction
c=inputStruct.numPeriods(3); % in z direction 

Ns=inputStruct.Ns; % number of sampling points
isocap=inputStruct.isocap;
levelset=inputStruct.levelset;
typeStr=inputStruct.surfaceCase;


%% Generation of points

%Create coordinates
rRange=linspace(0,(r2-r1),Ns);
phiRange=linspace(0,theta,Ns);
zRange=linspace(0,z,Ns);
[r_aux,phi,Z]=meshgrid(rRange,phiRange,zRange);

% Compute r_aux to actual radii
r= (r2)*r_aux + (r1)*(1- r_aux);

r=(r./r2);
% Z=((Z./abs(kz)).*c);

% compute Cartesian coordinates for grid points
Y = r .* cos(phi);
X = r .* sin(phi);

% Scaling coordinates
kx = (2*pi*a)*r_aux; 
ky = b*phi;
kz = (2*pi*c)*Z;

%% Generation of level sets and exporting
%Calculate 3D image data

%Evaluate metric on coordinates
switch typeStr
    case 'p' %Schwarz P
        S=(cos(kx)+cos(ky)+cos(kz));
    case 'd' %Schwarz D
        S=(sin(kx).*sin(ky).*sin(kz))...
            +(sin(kx).*cos(ky).*cos(kz))...
            +(cos(kx).*sin(ky).*cos(kz))...
            +(cos(kx).*cos(ky).*sin(kz));
    case 'g' %Schoen Gyroid
        S=(sin(kx).*cos(ky))+(sin(ky).*cos(kz))+(cos(kx).*sin(kz));
    case 'n' %Neovius
        S=3*(cos(kx)+ cos(ky)+ cos(kz))+ (4*cos(kx).*cos(ky).*cos(kz));
    case 'w'
        S=2*(cos(kx).*cos(ky)+cos(kz).*cos(kx)+cos(ky).*cos(kz))-(cos(2*kx)+cos(2*ky)+cos(2*kz));
    case 'pw'
        S=(4.*(cos(kx).*cos(ky)+cos(ky).*cos(kz)...
            +cos(kz).*cos(kx))-3.*cos(kx).*cos(ky).*cos(kz))+2.4;
    otherwise
        error('unknown surface type requested')
end

S=reshape(S,size(X));

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
% 
% % Visualizing geometry
% cFigure; hold on;
% gpatch(f,v,c,'none', 1);
% axisGeom(gca,20); axis on;
% colormap gjet; icolorbar;
% camlight Right;
% drawnow;


%% Collect output
varargout{1}=f;
varargout{2}=v;
varargout{3}=c;
varargout{4}=S;
varargout{5}=X;
varargout{6}=Y;
varargout{7}=Z;
varargout{8}=r;
varargout{9}=phi;

end
