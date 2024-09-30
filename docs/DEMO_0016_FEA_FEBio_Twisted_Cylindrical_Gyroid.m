%% DEMO_0015_FEA_ABAQUS_Twisted_Cylindrical_Gyroid
% This is a demo for:
% 
% * Building a cylinder containing a radial Gyroid lattice structure
% * Defining boundary conditions for finite element analysis
% * Automatically generating an ABAQUS .inp file
% * Triggering finite element analysis using ABAQUS
% * Importing and visualising the results 
%%

clear; close all; clc;

%%
%
%  Change log:
% 2023 VM created
% 2024/09/03 KMM Added header description
% 2024/09/24 KMM Edited to enable FEBio based simulation
%% Plot settings

cMap=spectral(250);
faceAlpha1=1;
faceAlpha2=0.65;
edgeColor1='none';
edgeColor2='none';
fontSize=25; 
pColors=gjet(6);
markerSize=20; 

%% Control Parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'assets','temp');
if ~exist(savePath,'dir')
    mkdir(savePath)
end

% Defining file names
febioFebFileNamePart='twistedCylinder_FEA';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress sigma_z

% Define applied displacement & direction
displacementMagnitude = -0.2; %

% Material parameter set
E_youngs=5000;
v_poisson=0.33;

% Set parameters for individual gyroid
% NOTE: Not all geometry parameters yield a valid mesh. If the image resolution 
% is poor, or if the point spacing is too large, the mesh may be invalid. 

radiusInner = 0.3; 
radiusOuter = 1.2; 
height = 2; 
Dim=[radiusInner, radiusOuter, 2*pi, height]; % size vector
Ns=150; % number of sampling points
numPeriods=[2 12 2]; %Number of periods in each direction
levelset=0.5 ; %Isosurface level

pointSpacing=(radiusOuter-radiusInner)/40;
tolDir=pointSpacing/5; %Tolerance for detecting sides after remeshing

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

runMode='internal';% 'internal' or 'external'

%% Compute gyroid sample in cilyndrical arrangment

% Set parameters for individual gyroid
inputStruct.size= Dim; % characteristic length
inputStruct.Ns=Ns; % number of sampling points
inputStruct.isocap=1; %Option to cap the isosurface
inputStruct.surfaceCase='g'; %Surface type
inputStruct.numPeriods=numPeriods; %Number of periods in each direction
inputStruct.levelset=levelset ; %Isosurface level

%% Evaluating cylindrical gyroid function 

[~,~,~,S,X,Y,Z]=TPMS_LCS(inputStruct);

%% Apply deformation matrix on grids 

bendAngle = pi; %twist angle
a=linspace(0,bendAngle,size(Z,3));

Xp=X; Yp=Y; Zp=Z;
for q=1:1:size(Z,3)
    R = euler2DCM([0 0 a(q)]);
    x = X(:,:,q);
    y = Y(:,:,q);
    z = Z(:,:,q);
    v = [x(:) y(:) z(:)];
    %     vp = v*R;
    vp = (R*v')';

    xp = reshape(vp(:,1),size(x));
    yp = reshape(vp(:,2),size(x));
    zp = reshape(vp(:,3),size(x));

    Xp(:,:,q)=xp;
    Yp(:,:,q)=yp;
    Zp(:,:,q)=zp;
end

%% Create new twisted grids 

V=[X(:) Y(:) Z(:)];  %original grids
Vp=[Xp(:) Yp(:) Zp(:)]; %deformed grids

% Visualization
h1=cFigure; hold on; 
title ('Original and Deformed Vertices', 'FontSize',fontSize);
plotV(V,'ko','MarkerSize',5);
plotV(Vp,'r*','MarkerSize',3);
legend({'Origional Grid','Deformed Grid'}, 'FontSize',20);
axisGeom; 
drawnow; 

%% Construct iso-surface

Sn=S;
[Fi,Vi] = isosurface(Xp,Yp,Zp,Sn,levelset);
[Fc,Vc] = isocaps(Xp,Yp,Zp,Sn,levelset);

[F,V] = joinElementSets({Fi,Fc},{Vi,Vc});
[F,V] = mergeVertices(F,V);
F=fliplr(F);

% Visualizing geometry
cFigure; hold on;
% title('GY-Scale-2','FontSize',fontSize);
gpatch(F,V,[0.75 0.75 0],'none', 1);
axisGeom(gca,fontSize); axis off;
colormap gjet; icolorbar; 
camlight left;
drawnow;

%% Exporting STL file 

% stlStruct.solidNames={'a'}; %names of parts
% stlStruct.solidFaces={F}; %Faces
% stlStruct.solidVertices={V}; %Vertices
% stlStruct.solidNormals={[]}; %Face normals (optional)
% 
% %Set main folder and fileName
% defaultFolder=fileparts(fileparts(mfilename('fullpath')));
% pathName=fullfile(defaultFolder, 'assets', 'STL');
% fileName=fullfile(pathName,'cylindrical_twisted_gyroid.stl');
% 
% export_STL_txt(fileName,stlStruct);

%% Remesh using geomgram

optionStruct.pointSpacing=pointSpacing;
optionStruct.pre.max_hole_area=0; %Max hole area for pre-processing step
optionStruct.pre.max_hole_edges=0; %Max number of hole edges for pre-processing step
[F,V]=ggremesh(F,V,optionStruct);
C=zeros(size(F,1),1);

% Visualizing geometry
cFigure; hold on;
title('Geogram remeshed','FontSize',fontSize);
gpatch(F,V,'w','k',1);
axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Get region point for meshing

V_region = getInnerPoint(F,V);

% Visualizing geometry
cFigure; hold on;
title('Geogram remeshed','FontSize',fontSize);
gpatch(F,V,'w','none',0.5);
plotV(V_region,'r.','MarkerSize',markerSize)
axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Tetrahedral meshing using tetgen (see also |runTetGen|)
% Create tetgen input structure
inputStruct.stringOpt='-pq1.2AaY';
inputStruct.Faces=F;
inputStruct.Nodes=V;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=C; %Face boundary markers
inputStruct.regionPoints=V_region; %region points
inputStruct.regionA=2*tetVolMeanEst(F,V);
inputStruct.minRegionMarker=2; %Minimum region marker

% Mesh model using tetrahedral elements using tetGen
[meshOutput]=runTetGen(inputStruct); %Run tetGen

% Access model element and patch data
Fb=meshOutput.facesBoundary;
F=meshOutput.faces;
Cb=meshOutput.boundaryMarker;
V=meshOutput.nodes;
CE=meshOutput.elementMaterialID;
E=meshOutput.elements;

%% Visualizing mesh using |meshView|, see also |anim8|

meshView(meshOutput);

%% Defining node labels

C_vertex=zeros(size(V,1),1);

logic1=V(:,3)>(Dim(1,4)-tolDir);
logic2=V(:,3)<tolDir;

C_vertex(logic1)=max(C_vertex(:))+1;
C_vertex(logic2)=max(C_vertex(:))+1;

% Visualizing vertex/node labels
hf=cFigure;
title('Boundary Nodes','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,'w','none',1);
scatterV(V,10,C_vertex,'filled');

axisGeom(gca,fontSize);
colormap gjet; icolorbar;
camlight headlight;
drawnow;

%% Defining the boundary conditions

bcSupportList=find(C_vertex==2); %Bottom vertices
bcPrescribeList=find(C_vertex==1); %Top vertices

% Visualizing boundary conditions. 
% Markers plotted on the semi-transparent model denote the nodes in the
% various boundary condition lists.

hl=cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
gpatch(Fb,V,'kw','none',0.5);

hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize);
hl(2)=plotV(V(bcPrescribeList,:),'r*','MarkerSize',markerSize/3);
legend(hl,{'BC support','BC prescribe'});

axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='4.0'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

%Control section
febio_spec.Control.analysis='STATIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.qn_method.max_ups=max_ups;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;

%Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.E=E_youngs;
febio_spec.Material.material{1}.v=v_poisson;

% Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='tet4'; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix
 
% -> NodeSets
nodeSetName1='bcSupportList';
nodeSetName2='bcPrescribeList';

febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList);

febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcPrescribeList);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_xyz';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof=1;
febio_spec.Boundary.bc{1}.y_dof=1;
febio_spec.Boundary.bc{1}.z_dof=1;

febio_spec.Boundary.bc{2}.ATTR.name='zero_displacement_xy';
febio_spec.Boundary.bc{2}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{2}.x_dof=1;
febio_spec.Boundary.bc{2}.y_dof=1;
febio_spec.Boundary.bc{2}.z_dof=0;

febio_spec.Boundary.bc{3}.ATTR.name='prescibed_displacement_z';
febio_spec.Boundary.bc{3}.ATTR.type='prescribed displacement';
febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{3}.dof='z';
febio_spec.Boundary.bc{3}.value.ATTR.lc=1;
febio_spec.Boundary.bc{3}.value.VAL=displacementMagnitude;
febio_spec.Boundary.bc{3}.relative=0;

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.name='LC_1';
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{1}.extend='CONSTANT';
febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0 0; 1 1];

febio_spec.Output.plotfile.var{end+1}.ATTR.type='right stretch';

%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='s1;s2;s3';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

% Plotfile section
febio_spec.Output.plotfile.compression=0;

%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window. 

%%1
% |febView(febio_spec); %Viewing the febio file|1

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function. 

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
%system(['gedit ',febioFebFileName,' &']);

%% Running the FEBio analysis
% To run the analysis defined by the created FEBio input file the
% |runMonitorFEBio| function is used. The input for this function is a
% structure defining job settings e.g. the FEBio input file name. The
% optional output runFlag informs the user if the analysis was run
% succesfully. 

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.runMode=runMode;
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%%     
% Importing nodal displacements from a log file
dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),0,1);

%Access data
N_disp_mat=dataStruct.data; %Displacement
timeVec=dataStruct.time; %Time

%Create deformed coordinate set
V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);

%% 
% Plotting the simulated results using |anim8| to visualize and animate
% deformations 

DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude
    
% Create basic view and store graphics handle to initiate animation
hf=cFigure; %Open figure  
gtitle([febioFebFileNamePart,': Press play to animate']);
title('Displacement magnitude [mm]','Interpreter','Latex')
hp=gpatch(Fb,V_DEF(:,:,end),DN_magnitude,'k',1,0.5); %Add graphics object to animate
hp.FaceColor='interp';
gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object

axisGeom(gca,fontSize); 
colormap(cMap); colorbar;
caxis([0 max(DN_magnitude)]); caxis manual;   
axis(axisLim(V_DEF)); %Set axis limits statically    
view(140,30);
camlight headlight;        
    
% Set up animation features
animStruct.Time=timeVec; %The time vector    
for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
    DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); %Current displacement magnitude
            
    %Set entries in animation structure
    animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
    animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
    animStruct.Set{qt}={V_DEF(:,:,qt),DN_magnitude}; %Property values for to set in order to animate
end        
anim8(hf,animStruct); %Initiate animation feature    
drawnow;

%% 
% _*LatticeWorks footer text*_ 
% 
% License: <https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE>
% 
% Copyright (C) 2023 Mahtab Vafaeefar and the LatticeWorks contributors
