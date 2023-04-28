%--------------------------------- Newmark -----------------------------%
function [disp,vel,acc,M,C,K,fem] = FEM_Newmark(fem,E,V)
gama = fem.gamma; beta = fem.beta; % Neewmark-beta parameters
dt = fem.Tmax/fem.NStep+1; % Time increment
f = fem.FreeDofs;
a0 = 1/beta/dt^2;
a1 = gama/beta/dt;
a2 = 1/beta/dt;
a3 = 1/2/beta-1;
a4 = gama/beta-1;
a5 = dt/2*(gama/beta-2);
a6 = dt*(1-gama);
a7 = gama*dt;
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
Kn = K+a0*M+a1*C;
acc(f,1) = M(f,f)\(Ft(f,1)-K(f,f)*disp(f,1)-C(f,f)*vel(f,1));
for It = 2:fem.NStep+1
    Fn = Ft(:,It)+M*(a0*disp(:,It-1)+a2*vel(:,It-1)+a3*acc(:,It-1))...
        +C*(a1*disp(:,It-1)+a4*vel(:,It-1)+a5*acc(:,It-1));
    if It==2
        [disp(f,It),L,s] = SolveLinSys(Kn(f,f),Fn(f));
    else
        disp(f(s),It) = L'\(L\Fn(f(s)));
    end    
    acc(:,It) = a0*(disp(:,It)-disp(:,It-1))-a2*vel(:,It-1)-a3*acc(:,It-1);
    vel(:,It) = vel(:,It-1)+a6*acc(:,It-1)+a7*acc(:,It);
end
%% ------------------------------------- GLOBAL STIFFNESS AND MASS MATRICES
function [M,C,K] = GlobalMCK(fem,E,V)
K = sparse(fem.i,fem.j,E(fem.e).*fem.k0); % Assemble stiffness matrix
K = (K+K')/2; 
M = sparse(fem.i,fem.j,V(fem.e).*fem.m0); % Assemble mass matrix
M = (M+M')/2;
C = fem.Ar(1)*M + fem.Ar(2)*K;            % Compute damping matrix
%% ----------------------- SOLVE LINEAR SYSTEM USING CHOLESKY DECOMPOSITION
function [U, L, s] = SolveLinSys(K,F)
[L,~,s] = chol(K,'lower','vector');
U(s,:) =L'\(L\F(s,:));
%-------------------------------------------------------------------------%