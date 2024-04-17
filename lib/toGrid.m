function [rho_VG, X,Y,Z, DG_im, imOrigin] = toGrid(E,Fb,V,rho, voxelSize)

% [rho_VG, X,Y,Z, DG_im, imOrigin] = toGrid(E,Fb,V,rho, voxelSize)
%
% -----------------------------------------------------------------------
% This function is used to map an element (E) desnity distribution field (rho), to a
% desired voxelSize (voxelSize) grid. 
%
% _*Name*_ 
% 
% License: <hyperlink to license>
% Author: _Mahtab Vafaee_, <mahtab.vafaee@gmail.com>
%
%  Change log:
%  2024/03/01 MV Created
%  2024/04/17 MV Edited
% -----------------------------------------------------------------------
%%

%% an image on the grids

[M,G,~]=patch2Im(Fb,V,[],voxelSize, [],[],2); % converting to image
L= M==1; % selecting the inside voxels

imOrigin=G.origin;
[J,I,K]=meshgrid(1:1:size(M,2),1:1:size(M,1),1:1:size(M,3));
[X,Y,Z]=im2cart(I,J,K,voxelSize);

VG=[X(:) Y(:) Z(:)]; % grid vertices
V=V-imOrigin(ones(size(V,1),1),:); % shift the vertices

%% Distance image from boundary faces 

Vsm=patchCentre(Fb,V);% find the center of each face

[DG,indClosest]=minDist(VG,Vsm); % finding the closest Vsm to each grid
DG_im=reshape(DG,size(L));

Et=[Fb(indClosest,:) (1:numel(indClosest))'+size(V,1)]; % tetrahedral elements with on face from Fb & a vertix from VG
Vt=[V; VG]; % merging the verices

[VE,logicPositive]=tetVol(Et,Vt,0); % positive and negative volume of generated tet elements

logicPositive=reshape(logicPositive,size(L)); % logicPositive are outside

DG_im(~logicPositive)=-DG_im(~logicPositive); % apply negative image values to inner voxels

%% Filter elements inside the boundaries

elementIndexFound = 1:1:size(E,1); 
elementIndexFound (logicPositive)= NaN; % Exclude elements outside
logicInside = ~logicPositive;

%% Interpolation of rho onto grid

interpolationMethod = 2;
switch interpolationMethod
    case 1 % Neirest
        rho_VG=nan(size(VG,1),1);
        rho_VG(logicInside)=rho(elementIndexFound(logicInside));
        rho_VG=reshape(rho_VG, size(X));
    case 2 % Linear
        Vm=patchCentre(E,V); % center of the elements
        rho_VG=nan(size(VG,1),1);
        rho_VG(logicInside) = griddata(Vm(:,1),Vm(:,2),Vm(:,3),rho,X(logicInside),Y(logicInside),Z(logicInside),'linear'); %Interpolate on grid
    case 3 % Element nodal average tri-linear
        rho_V = faceToVertexMeasure(E,V,rho); % Average from elements to nodes

        % Shape function (=barycentric coordinate) based within element
        % tri-linear interpolation
        [E,V,~]=hex2tet(E,V,rho,4); % if the elements are hexahedron
        TR = triangulation(E,V); %Convert to "triangulation" type
        
        %Default options
        inputStruct.distCropOpt=1;
        inputStruct.chullCropOpt=1;
        inputStruct.waitbarOpt=1;
        inputStruct.toleranceMagnitude=1e-3; %
        [elementIndexFound,baryCentricCoordinate]=pointLocationTR(TR,VG,inputStruct); % Compute
        logicInside=~isnan(elementIndexFound); %if nan then the voxels are outside shape


        rho_VG=nan(size(VG,1),1);
        for indPoint = find(logicInside)'
            indElement =  elementIndexFound(indPoint);
            indNodes = E(indElement,:);
            rho_VG(indPoint) = sum(baryCentricCoordinate(indPoint,:) .* rho_V(indNodes)');
        end
    case 4 % Natural
        Vm=patchCentre(E,V); % center of the elements
        interpFunc_rho = scatteredInterpolant(Vm,rho,'natural','none'); %Create interpolator
        rho_VG=nan(size(VG,1),1);
        rho_VG(logicInside) = interpFunc_rho(VG(logicInside,:));
end
rho_VG=reshape(rho_VG, size(X));

end