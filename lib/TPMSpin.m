function [Fg,Vg,S]=TPMSpin(X,Y,Z,inputStruct)

% [Fg,Vg,S]=TPMSpin(X,Y,Z,inputStruct)
%
% -----------------------------------------------------------------------
% This function is used to TPMS or Spinodoid field values on input grids 
% (X,Y,Z), and retures surface S, faces (Fg), and vertices (Vg).
%
% _*Name*_ 
% 
% License: <hyperlink to license>
% Author: _Mahtab Vafaee_, <mahtab.vafaee@gmail.com>
%
%  Change log:
%  2024/02/06 MV Created
% -----------------------------------------------------------------------
%%
levelset=inputStruct.levelset;
freq=inputStruct.numPeriods;
latType=inputStruct.surfaceType;

% 
X=freq(1,1).*X;
Y=freq(1,2).*Y;
Z=freq(1,3).*Z;

switch latType
    case 'p'
        S=(cos(X)+cos(Y)+cos(Z));

    case 'g'
        S=(sin(X).*cos(Y))+(sin(Y).*cos(Z))+(cos(X).*sin(Z));

    case 'd'
        S=(sin(X).*sin(Y).*sin(Z))...
            +(sin(X).*cos(Y).*cos(Z))...
            +(cos(X).*sin(Y).*cos(Z))...
            +(cos(X).*cos(Y).*sin(Z));

    case 'n'
        S=3*(cos(X)+ cos(Y)+ cos(Z))+ (4*cos(X).*cos(Y).*cos(Z));

    case 'w'
        S=2*(cos(X).*cos(Y)+cos(Z).*cos(X)+cos(Y).*cos(Z))-...
            (cos(2*X)+cos(2*Y)+cos(2*Z));

    case 'pw'
        S=(4.*(cos(X).*cos(Y)+cos(Y).*cos(Z)...
            +cos(Z).*cos(X))-3.*cos(X).*cos(Y).*cos(Z))+2.4;

    case 'spin'
        R=eye(3);
        thetas=freq;
        relativeDensity=inputStruct.relativeDensity;
        waveNumber=inputStruct.waveNumber;
        numWaves=1000;

        % Generate wave phase angles for GRF
        wavePhases = rand_angle([numWaves,1]); %2*pi*rand(numWaves,1);

        % Generate wave directions for GRF
        %array of all wave directions
        waveDirections = zeros(numWaves,3);
        for i=1:numWaves
            flag = true; %keep trying until candidate wave vector is found
            while(flag)
                % generate isotropically sampled candidate wave
                candidate = randn(1,3);
                candidate = candidate/norm(candidate);
                % check for allowed wave vector directions
                % angle along first axis
                angle1 = min(...
                    acosd(dot(candidate,R(1,:))),...
                    acosd(dot(candidate,-R(1,:))));
                % angle along second axis
                angle2 = min(...
                    acosd(dot(candidate,R(2,:))),...
                    acosd(dot(candidate,-R(2,:))));
                % angle along third axis
                angle3 = min(...
                    acosd(dot(candidate,R(3,:))),...
                    acosd(dot(candidate,-R(3,:))));
                % check
                if(any([angle1,angle2,angle3]<thetas))
                    flag = false;
                    break;
                end
            end
            waveDirections(i,:)=candidate;
        end
        % Evaluate GRF on sampling points
        S = zeros(size(X));
        for i=1:numWaves
            dotProduct = waveDirections(i,1)*X ...
                + waveDirections(i,2)*Y ...
                + waveDirections(i,3)*Z;

            S = S+sqrt(2/numWaves)*cos(dotProduct*waveNumber+wavePhases(i));
        end

        % Apply levelset
        levelset= sqrt(2)*erfinv(2*relativeDensity-1);
end

% Scaling coordinates
X=X./freq(1,1);
Y=Y./freq(1,2);
Z=Z./freq(1,3);

%Construct iso-surface
[Fg,Vg] = isosurface(X,Y,Z,S,levelset);

end