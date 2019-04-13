%% Prob 2 (perm.A)
clearvars
close all

%Perm.A
A=1;

q=[0.4,1.0];  %interp.1D, using all nodes at the upper boundary.
r=[0.4,0.3];  %interp.2D, using quadrilaterls.

f=@(x,y) exp(-1./(A+x.^2+y.^2));

eval('mesh8x8Trapz');

numNod=size(nodes,1);
numElem=size(elem,1); 

% Plot mesh and points q and p
plotQuadMesh(nodes,elem)
hold on
plot(q(1,1),q(1,2),'o',...
    'Marker','o','markerFaceColor','red','markerSize',6,...
    'color','black','lineWidth',1.25)
plot(r(1,1),r(1,2),'o',...
    'Marker','o','markerFaceColor','yellow','markerSize',6,...
    'color','black','lineWidth',1.25)

%Indexs of Nodes at the top boundary
idxTop=find(nodes(:,2) > 0.9);
idxTop=idxTop';

plot(nodes(idxTop,1),nodes(idxTop,2),'o',...
    'Marker','o','markerFaceColor','blue','markerSize',6,...
    'color','black','lineWidth',1.25)

tempNodTop = f(nodes(idxTop,1),nodes(idxTop,2));

%1D interpolation at point q
n=size(idxTop,2);
interpPol=polyfit(nodes(idxTop,1),tempNodTop,n-1);
interpQ=polyval(interpPol,q(1,1));

TempQ=f(q(1,1),q(1,2));
errInterpQ=abs(interpQ-TempQ);

%2D interpolation at point r
for e=1:numElem
    v1=nodes(elem(e,1),:);
    v2=nodes(elem(e,2),:);
    v3=nodes(elem(e,3),:);
    v4=nodes(elem(e,4),:);
    vertexs=[v1;v2;v3;v4];
    [alphas,isInside]=baryCoordQuad(vertexs,r);
    if (isInside >= 1)
        tempNodesElem=f(vertexs(:,1),vertexs(:,2));
        break;
    end 
end

interpR=alphas*tempNodesElem;
tempR=f(r(1,1),r(1,2));
errInterpR=abs(interpR-tempR);

%Output results
fprintf("%43s\n",'Error')
fprintf("Interp.Temp.at point R: %10.4e%12.4e\n",interpR,errInterpR)
fprintf("Interp.Temp.at point Q: %10.4e%12.4e\n",interpQ,errInterpQ)