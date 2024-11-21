function [graded_S]=multiMorph (X,Y,Z,Data)

% function [graded_S]=multiMorph (X,Y,Z,Data);
%
% -----------------------------------------------------------------------
% This function is used to compute multi-morphology surface from individual
% lattices.
%
% for |multiMorph|, (X,Y,Z) are the coordinates of grids 
% and  (Data) is the input structre with all individual surfaces and the
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

if ~isfield(Data,'T') || isempty (Data.T)
    Data.T=1; % Default transition type
end
transType=Data.T; % Type {(1=square_distance), (2=)}

% kappa controls the lengthscale of transition between lattices
% Higher kappa => faster transition
% Lower kappa => slower transition
kappa = Data.K;

l=length(Data.C);

% Initializing parameters size
sum_weights = zeros(size(Data.S{1}));
weights = cell(l, 1);
S = cell(l, 1);
center = cell(l, 1);
levelset = zeros(l, 1);

% Loop over the individual morphologies/lattices
for i= 1:l
    S{i}= Data.S{i};
    center{i} = Data.C{i};
    levelset(i) = Data.L{i};
    
    % Using Gaussian (a.k.a. radial basis functions) interpolation.
    % Computing the weights for each lattice evaluated on all grid points.

    if transType==1
        weights{i} = exp(-kappa * Squared_distance_from_point(X,Y,Z,center{i}));
    elseif transType==2
        weights{i} = exp(-kappa * ((X-center{i}(1,1)).^2 + (Y-center{i}(1,2)).^2));
    elseif transType==3 
        weights{i} = exp(-kappa * (((X.^2+Y.^2)-(center{i}(1,1)).^2 + center{i}(1,2)).^2));
    end

    % sum of weitghs must be 1
    sum_weights = sum_weights + weights {i};
end

graded_S = zeros (size(S{1}));
for i= 1:l
    weights{i} = weights {i}./sum_weights; % Normalizing the weights
    weighted_S = weights{i} .* (S{i} - levelset(i));

    graded_S = graded_S + weighted_S; % Interpolating using the above weights
end

end

%%
% _*LatticeWorks footer text*_
%
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
%
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors

