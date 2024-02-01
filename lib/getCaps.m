function [fc,vc,cc]=getCaps(X,Y,Z,S,levelset)

% function [fc,vc,cc]=getCaps(X,Y,Z,S,levelset);
%
% -----------------------------------------------------------------------
% This function is used to cap the boundaries of a generated surface.
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
    [fc,vc] = isocaps(X,Y,Z,S,levelset);     %Compute isocaps
    
    nc=patchNormal(fc,vc);
    cc=zeros(size(fc,1),1);
    cc(nc(:,1)<-0.5)=1;
    cc(nc(:,1)>0.5)=2;
    cc(nc(:,2)<-0.5)=3;
    cc(nc(:,2)>0.5)=4;
    cc(nc(:,3)<-0.5)=5;
    cc(nc(:,3)>0.5)=6;    
end