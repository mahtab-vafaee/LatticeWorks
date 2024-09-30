function [] = Hybrid_vizualize(varargin)

% function [h]=Hybrid_vizualize(F,V,C,map,center_V);
%
% -----------------------------------------------------------------------
% This function is a customized version of |gpatch| function, specifically
% designed for multi-morphology and hybrid lattices visualization.
%
% for |Hybrid_vizualize| are the faces (F), the vertices (V), the color
% description (C), the colormap description (map), the vertices of centre
% locations (center_V).
%
% _*LatticeWorks*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
%
% Author: _Mahtab Vafaee_, <mahtab.vafaee@gmail.com>
%
%  Change log:
%  2024/1/31 MV Created
%  2024/2/1 MV Fixed the variable input cases
% -----------------------------------------------------------------------
%% Parse inputs

switch nargin
    case 3
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
    case 5
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        map=varargin{4}; %Plotting map data
        center_V=varargin{5}; %Vertices of location centers
    otherwise
        error('Wrong number of inputs');
end

%% Using grouping to keep only largest group
groupOptStruct.outputType='label';
[G,~,groupSize]=tesgroup(F,groupOptStruct); %Group connected faces
[~,indKeep]=max(groupSize); %Index of largest group

%Keep only largest group
F=F(G==indKeep,:); %Trim faces
C=C(G==indKeep,:); %Trim color data
[F,V]=patchCleanUnused(F,V); %Remove unused nodes

%% Visualize surface
cFigure;
switch nargin
    case 3
        C=ones(size(F,1),1)*[0.75 0.75 0];%face colors
        map=[0.75 0.75 0];
        gpatch(F,V,C,'none', 1); %visualize faces

    case 5
        if isempty (map)
            map=ones(7,1)*[0.75 0.75 0];
        end
        if isempty (center_V)
            gpatch(F,V,C,'none', 1); %visualize faces
        end

        gpatch(F,V,C,'none', 0.5); %visualize surface
        hold on;
        %Center vertices
        for i=1:size(center_V,1)
            plotV(center_V(i,:),'r*', 'LineWidth',4.5,'markerSize',30);
            hold on;
        end
end

%% Visualization setting

axisGeom; camlight headlight;
colormap (map); icolorbar;
grid off; axis on;
gdrawnow;

end
%% 
% _*LatticeWorks footer text*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
% 
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
