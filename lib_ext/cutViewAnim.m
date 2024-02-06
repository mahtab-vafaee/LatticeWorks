function cutViewAnim(F,V)
%% function cutViewAnim(F,V), from cutViewAnim8.m GIBBON

FV=patchCentre(F,V)/5;
C=sin(FV(:,1))+sin(FV(:,2))+sin(FV(:,3));
snapTolerance=mean(patchEdgeLengths(F,V))/50; %Snapping tolerance
n=vecnormalize([0 0 1]); %Plane normal vector

P=mean(V,1); %Point on plane

[Fc,Vc,Cc,logicSide]=triSurfSlice(F,V,C,P,n,snapTolerance);
Eb=patchBoundary(Fc(logicSide,:));

hf=cFigure; 

hold on; 
hp1=gpatch(Fc(logicSide,:),Vc,'w','none',1,3);
hp2=gpatch(Fc(~logicSide,:),Vc,'w','none',0.25);
hp3=gpatch(Eb,Vc,'none','b',1,4);

axisGeom; axis manual; camlight headligth;
view(-90,0); zoom(1.25);axis off; 
gdrawnow; 

nSteps=75; %Number of animation steps
animStruct.Time=linspace(0,1,nSteps); %Create the time vector
z=linspace(min(V(:,3)),max(V(:,3)),nSteps);
for q=1:1:nSteps    
    P=[0 0 z(q)];    
    [Fc,Vc,~,logicSide,Eb]=triSurfSlice(F,V,C,P,n,snapTolerance);

    %Set entries in animation structure
    animStruct.Handles{q}=[hp1 hp1 hp2 hp2 hp3 hp3]; %Handles of objects to animate
    animStruct.Props{q}={'Vertices','Faces',...
                         'Vertices','Faces',...   
                         'Vertices','Faces',... 
                         'Vertices','Faces','CData'... 
                         }; %Properties of objects to animate
    animStruct.Set{q}={Vc,Fc(logicSide,:),Vc,Fc(~logicSide,:),Vc,Eb,Vc,Fc,double(logicSide)}; %Property values for to set in order to animate
end
anim8(hf,animStruct);

end