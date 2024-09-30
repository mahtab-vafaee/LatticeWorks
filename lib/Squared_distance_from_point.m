function [sq_dist] = Squared_distance_from_point(X,Y,Z,point)

% function [sq_dist]=Squared_distance_from_point(X,Y,Z,point);
%
% -----------------------------------------------------------------------
% This function is used to compute squared distance from a point
%
% for |Squared_distance_from_point|, (X,Y,Z) are the coordinates of grids 
% and  (point) is the requested point to compute distance from.
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
%
% X, Y, Z are 3D matrices
% point is a 1x3 vector
sq_dist = (X-point(1)).^2 + (Y-point(2)).^2 + (Z-point(3)).^2;

end
%% 
% _*LatticeWorks footer text*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
% 
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
