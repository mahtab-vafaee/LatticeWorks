function [f,v,c]=getSurface(X,Y,Z,S,levelset)

% function [f,v,c]=getSurface(X,Y,Z,S,levelset);
%
% -----------------------------------------------------------------------
% This function is used to draw an isosurface over (X,Y,Z) grid with field
% variable (S), at specified levelset (levelset).
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
    [f,v] = isosurface(X,Y,Z,S,levelset);
    c=zeros(size(f,1),1);
end

