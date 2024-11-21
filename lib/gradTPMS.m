function [S,X,Y,Z] = gradTPMS(inputStruct)
%
% function [S,X,Y,Z] = gradTPMS(inputStruct);
%
% -----------------------------------------------------------------------
% This function is used to compute graded TPMS with cell size or 
% volume fraction gradient.
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

%Get parameters from input structure
L = inputStruct.L; % characteristic length
Ns = inputStruct.Ns; % number of sampling points
k = inputStruct.numPeriods;
gradType = inputStruct.gradType;
surfaceCase = inputStruct.surfaceCase;


%Create coordinates
xMin=0; xMax= 2*pi*k(1,1);
yMin=0; yMax= 2*pi*k(1,2);
zMin=0; zMax= 2*pi*k(1,3);

xRange=linspace(xMin,xMax,Ns);
yRange=linspace(yMin,yMax,Ns);
zRange=linspace(zMin,zMax,Ns);
[X,Y,Z]=meshgrid(xRange,yRange,zRange);

switch gradType
    case 'levelSet' % volume fraction gradient, Figure 1 (a,b)
        %Calculate 3D image data

        minL = inputStruct.GF(1,1);
        maxL = inputStruct.GF(1,2);

        switch surfaceCase
            case 'p' %Schwarz P
                S=(cos(X)+cos(Y)+cos(Z));
            case 'd' %Schwarz D
                S=(sin(X).*sin(Y).*sin(Z))...
                    +(sin(X).*cos(Y).*cos(Z))...
                    +(cos(X).*sin(Y).*cos(Z))...
                    +(cos(X).*cos(Y).*sin(Z));
            case 'g' %Schoen Gyroid
                S=(sin(X).*cos(Y))+(sin(Y).*cos(Z))+(cos(X).*sin(Z));
            case 'n' %Neovius
                S=3*(cos(X)+ cos(Y)+ cos(Z))+ (4*cos(X).*cos(Y).*cos(Z));
            case 'w'
                S=2*(cos(X).*cos(Y)+cos(Z).*cos(X)+cos(Y).*cos(Z))-(cos(2*X)+cos(2*Y)+cos(2*Z));
            case 'pw'
                S=(4.*(cos(X).*cos(Y)+cos(Y).*cos(Z)...
                    +cos(Z).*cos(X))-3.*cos(X).*cos(Y).*cos(Z))+2.4;
            otherwise
                error('unknown surface type requested')
        end

        S=reshape(S,size(X));

        % levelset gradient
        GF=X; % Use x-dir for now
        GF=GF-min(GF(:)); % 0-...
        GF=GF./max(GF(:)); % 0-1

        GF=GF*((-1/maxL)-(-1/minL)); % 0-2.5
        GF=GF + (-1/minL); % 0.8333-3.3333

        S=S.*GF;


    case 'cellSize' % cell size gradient, Figure 1 (d,e)
        % Calculate gradient frequency
        m = inputStruct.GF;
        phaseShift = inputStruct.phaseShift;

        K1= (m-1)/(xMax-xMin);
        C1= (xMin*K1)+1;
        C0= 0.5*K1*(xMin)^2;

        a = K1/2*X+C1+C0/X;
        b = K1*X+C1;
        c = K1*X+C1;

        %Calculate 3D image data

        switch surfaceCase
            case 'p' %Schwarz P
                S=(cos(a.*(X-phaseShift))+cos(b.*(Y-phaseShift))+cos(c.*(Z-phaseShift)));
            case 'd' %Schwarz D
                S=(sin(a.*(X-phaseShift)).*sin(b.*(Y-phaseShift)).*sin(c.*(Z-phaseShift)))...
                    +(sin(a.*(X-phaseShift)).*cos(b.*(Y-phaseShift)).*cos(c.*(Z-phaseShift)))...
                    +(cos(a.*(X-phaseShift)).*sin(b.*(Y-phaseShift)).*cos(c.*(Z-phaseShift)))...
                    +(cos(a.*(X-phaseShift)).*cos(b.*(Y-phaseShift)).*sin(c.*(Z-phaseShift)));
            case 'g' %Schoen Gyroid
                S=( sin(a.*(X-phaseShift)).*cos(b.*(Y-phaseShift)))+(sin(b.*(Y-phaseShift)).*...
                    cos(c.*(Z-phaseShift)))+(cos(a.*(X-phaseShift)).*sin(c.*(Z-phaseShift)) );
            case 'n' %Neovius
                S=3*(cos(a.*(X-phaseShift))+ cos(b.*(Y-phaseShift))+cos(c.*(Z-phaseShift)))...
                    + (4*cos(a.*(X-phaseShift)).*cos(b.*(Y-phaseShift)).*cos(c.*(Z-phaseShift)));
            case 'w'
                S=2*(cos(a.*(X-phaseShift)).*cos(Y)+cos(c.*(Z-phaseShift)).*cos(a.*(X-phaseShift))...
                    +cos(b.*(Y-phaseShift)).*cos(c.*(Z-phaseShift)))-(cos(2*a.*(X-phaseShift))+...
                    cos(2*b.*(Y-phaseShift))+cos(2*c.*(Z-phaseShift)));
            case 'pw'
                S=(4.*(cos(a.*(X-phaseShift)).*cos(b.*(Y-phaseShift))...
                    +cos(b.*(Y-phaseShift)).*cos(c.*(Z-phaseShift))...
                    +cos(c.*(Z-phaseShift)).*cos(a.*(X-phaseShift)))...
                    -3.*cos(a.*(X-phaseShift)).*cos(b.*(Y-phaseShift)).*cos(c.*(Z-phaseShift)))+2.4;
            otherwise
                error('unknown surface type requested')
        end

        S=reshape(S,size(X));
end

%% Scaling coordinates
switch length(L)
    case 1
        X=((X./abs(xMax-xMin)).*L);
        Y=((Y./abs(yMax-yMin)).*L);
        Z=((Z./abs(zMax-zMin)).*L);

    case 3
        X=((X./max(X(:))).*L(1,1));
        Y=((Y./max(Y(:))).*L(1,2));
        Z=((Z./max(Z(:))).*L(1,3));
end

end

%%
% _*LatticeWorks footer text*_
%
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
%
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
