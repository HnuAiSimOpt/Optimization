%--------------------------------- MSM -----------------------------%
function [disp,vel,acc,M,C,K,fem] = FEM_MSM(fem,E,V,Nm)
fr = fem.FreeDofs;
n = 3*fem.NNode;
disp = zeros(3*fem.NNode,fem.NStep+1); vel = disp; acc = disp;
disp(:,1) = fem.u0; vel(:,1) = fem.v0; % Initial conditions
[M,C,K] = GlobalMCK(fem,E,V); % Compute stiffness and mass matrices
if ~isempty(fem.Mass)
  DOF = [3*fem.Mass(:,1)-2;3*fem.Mass(:,1)-1;3*fem.Mass(:,1)];
  M(DOF,DOF)= M(DOF,DOF) + [fem.Mass(:,2);fem.Mass(:,2)];
end
if ~isempty(fem.ag)
  fem.Fa = -M*ones(3*fem.NNode,1).*fem.ag;
  Ft = fem.Fext+fem.Fa;
else
  Ft = fem.Fext;
end
Omiga0 = 2*pi;
[Vec,D] = eigs(K(fr,fr),M(fr,fr),Nm,'SM');
[~,nm] = sort(D);   Vec = Vec(:,nm);
V1 = zeros(length(fr),Nm); VN = V1;
for Im = 1:Nm
    V1(:,Im) = Vec(:,Im)/Vec(1,Im);
end
Mp = diag(diag(V1'*M(fr,fr)*V1));
for Im = 1:Nm
    VN(:,Im) = V1(:,Im)/sqrt(Mp(Im,Im)) ;
end
Omiga = diag(sqrt(VN'*K(fr,fr)*VN));
disY = VN'*Ft(fr,:)./(Omiga.^2-Omiga0^2);
for Im = 1:Nm
    disp(fr,:) = VN(:,Im)*disY(Im,:)- ...
        VN(:,Im)*(VN(:,Im)'*Ft(fr,:))/Omiga(Im)^2;
end
disp(fr,:) = SolveLinSys(K(fr,fr),Ft(fr,:)) + disp(fr,:);
%% ------------------------------------- GLOBAL STIFFNESS AND MASS MATRICES
function [M,C,K] = GlobalMCK(fem,E,V)
K = sparse(fem.i,fem.j,E(fem.e).*fem.k0); % Assemble stiffness matrix
K = (K+K')/2; 
M = sparse(fem.i,fem.j,V(fem.e).*fem.m0); % Assemble mass matrix
M = (M+M')/2;
C = fem.Ar(1)*M + fem.Ar(2)*K;            % Compute damping matrix
%% ----------------------- SOLVE LINEAR SYSTEM USING CHOLESKY DECOMPOSITION
function [U] = SolveLinSys(K,F)
[L,~,s] = chol(K,'lower','vector');
U(s,:) =L'\(L\F(s,:));
%-------------------------------------------------------------------------%