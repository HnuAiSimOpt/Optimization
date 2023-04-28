%--------------------------- PreComputations -----------------------------%
function fem = PreComputations(fem)
% Nodes ID of per element
V = cell(fem.NElem,1);
for el = 1:fem.NElem
    V{el} = unique([fem.Face{fem.Element{el}}]);
end
%Compute triplets used to assemble stiffness matrix and load vector
fem.ElemNDof = 3*cellfun(@length,V); % # of DOFs per element
fem.i = zeros(sum(fem.ElemNDof.^2),1);
fem.j=fem.i; fem.e=fem.i; fem.k0 = fem.i; fem.m0 = fem.i;
fem.Fa0 = zeros(sum(fem.ElemNDof),1); 
fem.eDof = fem.Fa0; %fem.DofE = fem.eDof; fem.Fa0=fem.eDof;
fem.iK0=fem.i; fem.jK0=fem.iK0;
indexK = 0; indexF = 0; IniI=0; IniJ=0;
for el = 1:fem.NElem     
  NDof = fem.ElemNDof(el); NKE = NDof^2; eNode = V{el};
  eDof = reshape([3*eNode-2;3*eNode-1;3*eNode],NDof,1); 
  %Triplets for stiffness matrix
  I = repmat(eDof ,1,NDof); J = I';
  fem.i(indexK+1:indexK+NKE) = I(:);
  fem.j(indexK+1:indexK+NKE) = J(:); 
  fem.e(indexK+1:indexK+NKE) = el;
  fem.DofE(indexF+1:indexF+NDof) = el; 
  fem.eDof(indexF+1:indexF+NDof) = eDof;
  %Local stiffness and mass matrices
  if ~fem.Reg || el==1,  [Ke, Me] = LocalKM3D(fem,fem.ShapeFnc{el},eNode); end
  fem.k0(indexK+1:indexK+NKE) = Ke(:);
  fem.m0(indexK+1:indexK+NKE) = Me(:);
  fem.Fa0(indexF+1:indexF+NDof) = Me*ones(NDof,1);
  %Element DOFs 
  I = repmat((1:NDof)',1,NDof); J=I';
  fem.iK0(indexK+1:indexK+NKE) = I(:)+IniI;
  fem.jK0(indexK+1:indexK+NKE) = J(:)+IniJ;
  IniI = IniI+NDof;
  IniJ = IniJ+NDof;
  indexK = indexK+NDof^2; indexF = indexF+NDof; 
end

fem.m0 = fem.rho.*fem.m0;
if ~isempty(fem.ag), fem.Fa0 = -fem.Fa0.*fem.ag;
else, fem.Fa0 = fem.Fa0.*zeros(1,fem.NStep+1); end
%Compute external load vector and initialize inertia force vector
NLoad = size(fem.Load,1);
fem.Fext = zeros(3*fem.NNode,fem.NStep+1);  %External load vector
fem.Fa = 0.*fem.Fext;                       %Inertia force vector
if NLoad>0
 fem.Fext(3*fem.Load(1:NLoad,1)-2,:) = fem.Load(1:NLoad,2:3:end);  %x-crdnt
 fem.Fext(3*fem.Load(1:NLoad,1)-1,:) = fem.Load(1:NLoad,3:3:end);  %y-crdnt
 fem.Fext(3*fem.Load(1:NLoad,1),:)   = fem.Load(1:NLoad,4:3:end);  %z-crdnt
end
%Obtain fixed DOFs and free DOFs
NSupp = size(fem.Supp,1);
FixedDofs = [fem.Supp(1:NSupp,2).*(3*fem.Supp(1:NSupp,1)-2);
             fem.Supp(1:NSupp,3).*(3*fem.Supp(1:NSupp,1)-1)
             fem.Supp(1:NSupp,4).*(3*fem.Supp(1:NSupp,1))];
FixedDofs = FixedDofs(FixedDofs>0);
AllDofs   = 1:3*fem.NNode;
fem.FreeDofs = setdiff(AllDofs,FixedDofs);
%-------------------------------------- ELEMENT STIFFNESS AND MASS MATRICES
function [Ke,Me] = LocalKM3D(vem,ShapeFnc,eNode)
nn = length(eNode);
D=vem.E0/((1+vem.Nu0)*(1-2*vem.Nu0))*[1-vem.Nu0 vem.Nu0 vem.Nu0 0 0 0;
                                      vem.Nu0 1-vem.Nu0 vem.Nu0 0 0 0;
                                      vem.Nu0 vem.Nu0 1-vem.Nu0 0 0 0;
                                                0 0 0 1/2-vem.Nu0 0 0;
                                                0 0 0 0 1/2-vem.Nu0 0;
                                                0 0 0 0 0 1/2-vem.Nu0];
Ke=zeros(3*nn,3*nn);Me = Ke;
W=ShapeFnc.W;
N = zeros(3,3*nn);
N(1,1:3:end) = ShapeFnc.N;
N(2,2:3:end) = N(1,1:3:end);
N(3,3:3:end) = N(1,1:3:end);
Proj=zeros(3*nn);
Proj(1:3:end,1:3:end)=ShapeFnc.P;
Proj(2:3:end,2:3:end)=ShapeFnc.P;
Proj(3:3:end,3:3:end)=ShapeFnc.P;   
B=zeros(6,3*nn); 
B(1,1:3:3*nn) = ShapeFnc.dNdX(:,1)';
B(2,2:3:3*nn) = ShapeFnc.dNdX(:,2)';
B(3,3:3:3*nn) = ShapeFnc.dNdX(:,3)';
B(4,1:3:3*nn) = ShapeFnc.dNdX(:,2)';
B(4,2:3:3*nn) = ShapeFnc.dNdX(:,1)';
B(5,2:3:3*nn) = ShapeFnc.dNdX(:,3)';
B(5,3:3:3*nn) = ShapeFnc.dNdX(:,2)';
B(6,1:3:3*nn) = ShapeFnc.dNdX(:,3)';
B(6,3:3:3*nn) = ShapeFnc.dNdX(:,1)';
Ke = Ke+B'*D*B*W;
Me = Me+N'*N*W;   %Mass matrix assuming unit density
alpha=vem.E0*(6-9*vem.Nu0)/9/(1-2*vem.Nu0)/(1+vem.Nu0);% VEM scaling factor
Ke=Ke+alpha*sum(W)^(1/3)*(eye(3*nn)-Proj')*(eye(3*nn)-Proj);
Ke=1/2*(Ke+Ke');
Me=Me+alpha*sum(W)^(1/3)*(eye(3*nn)-Proj')*(eye(3*nn)-Proj);
Me=1/2*(Me+Me');
