%% Prob.3 (perm.A)
clearvars
close all

%Perm.A
uA=3.0;
uB=8.0;
tol=0.1;
a1=3.0;
a0=2.0;
f=-3.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=2.0;
b=5.0;
L=b-a;

numElem=200;

%Geometry
h=L/numElem;           %elements of equal length!
nodes=(a:h:b)';
numNod=size(nodes,1);    
elem=zeros(numElem,2); %Connectivity matrix.
for e=1:numElem
    elem(e,:)=[e,e+1];
end
%Assembly of the global system.
K=zeros(numNod);
F=zeros(numNod,1);
u=zeros(numNod,1);
Fe=f*h/2.0*[1;1];
for e=1:numElem
    x1=nodes(elem(e,1),1); %1st. node of elem e.
    x2=nodes(elem(e,2),1); %2nd. node of elem e.
    Ke=a1*(x1*x1+x1*x2+x2*x2)/(3.0*h)*[1,-1;-1,1]+ ...
        a0*h/6.0*[2,1;1,2];
%    Ke=a1*(x2.^3-x1.^3)/(3.0*(h.^2))*[1,-1;-1,1]+ ...
%        a0*h/6.0*[2,1;1,2];
    rows=[elem(e,1),elem(e,2)];
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke;
    F(rows,1)=F(rows,1)+Fe;
end
fixedNod=[1,numNod];
freeNod=setdiff(1:numNod,fixedNod);

%Natural B.C.:
%set Q to zero. Remark: Q(1) and Q(end), .i.e., the
%components of Q corresponding to the fixed nodes are
%computed in the post-process (see below).
Q=zeros(numNod,1);

%Essential B.C.
u(1)=uA;
u(numNod)=uB;

%Reduced system.
Qm=Q(freeNod)-K(freeNod,fixedNod)*u(fixedNod);
Fm=F(freeNod)+Qm;
Km=K(freeNod,freeNod);
um=Km\Fm;
u(freeNod)=um;
    
%Post process.
Q=K*u-F;

%Mean value.
U=sum(u)/numNod; 
idx=find(abs(u-U)<tol);
N=size(idx,1);
fprintf('Solution:\n')

fprintf('(a) K(50,50) = %11.5e\n',K(50,50))
fprintf('(Hint K(20,21) = %11.5e)\n',K(20,21))
fprintf('(b) Mean val.of the approx.solution: ')
fprintf('<u> = %11.5e\n',U)
fprintf('(Hint: We know that u(25)= %11.5e)\n',u(25))
fprintf('(c) Number of nodes s.t. |u-<u>|< %.2f: N = %d\n',tol,N)