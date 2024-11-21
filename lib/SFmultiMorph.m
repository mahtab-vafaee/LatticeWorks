function [graded_S]=SFmultiMorph (inputStruct)
%
% function  [graded_S]=SFmultiMorph (inputStruct);
%
% -----------------------------------------------------------------------
% This function is used to compute multi-morphology surface using Sigmoid 
% function from individual lattices.
%
% for |SFmultiMorph|, (X,Y,Z) are the coordinates of grids 
% and  (inputStruct) is the input structre with all individual surfaces and the
% centers of regions to generate multi-Morphology lattice.
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

kappa = inputStruct.K;
G = inputStruct.G;
l=length(inputStruct.S);

% Initializing parameters size
S = cell(l, 1);
levelset = zeros(l, 1);

% Loop over the individual morphologies/lattices
for i= 1:l
    S{i}= inputStruct.S{i};
    levelset(i) = inputStruct.L{i};
end

%% Compute the weitgh functions

% Using Sigmoid Function interpolation.
% Computing the weights for each lattice evaluated on all grid points.
weights_A = 1/(1+exp(-kappa * G));
weights_B = (1-weights_A);
 
%% Interpolating using the above weights

graded_S =  weights_A .* (S{1} - levelset(1)) ...
    + weights_B.* (S{2} - levelset(2));

end 

%%
% _*LatticeWorks footer text*_
%
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
%
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors

