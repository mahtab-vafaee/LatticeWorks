function [f,v,c] = FV_arrange (f,v,fc,vc);

% function [f,v,c]=FV_arrange (f,v,fc,vc);
%
% -----------------------------------------------------------------------
% This function containes joining, merging, and cleaning functions for
% sorting faces and vertices after capping.
%
% for |FV_arrange|, (f,v) are the faces and vertices of the structure, and 
% (fc,vc) are faces and vertices from capping the structure.
%
% _*Name*_ 
% 
% License: <hyperlink to license>
% Author: _Mahtab Vafaee_, <mahtab.vafaee@gmail.com>
%
%  Change log:
%  2024/11/21 MV Created  
% -----------------------------------------------------------------------
%% 

c=zeros(size(f,1),1);

% Boilerplate code for preparing output for exporting/visualization
nc=patchNormal(fc,vc);
cc=zeros(size(fc,1),1);
cc(nc(:,1)<-0.5)=1;
cc(nc(:,1)>0.5)=2;
cc(nc(:,2)<-0.5)=3;
cc(nc(:,2)>0.5)=4;
cc(nc(:,3)<-0.5)=5;
cc(nc(:,3)>0.5)=6;

% Join sets
[f,v,c]=joinElementSets({f,fc},{v,vc},{c,cc});

% Merge nodes
[f,v]=mergeVertices(f,v);

% Check for unique faces
[~,indUni,~]=unique(sort(f,2),'rows');
f=f(indUni,:); %Keep unique faces
c=c(indUni);

% Remove collapsed faces
[f,logicKeep]=patchRemoveCollapsed(f);
c=c(logicKeep);

% Remove unused points
[f,v]=patchCleanUnused(f,v);

% Invert faces
f=fliplr(f);

end

%%
% _*LatticeWorks footer text*_
%
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
%
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
