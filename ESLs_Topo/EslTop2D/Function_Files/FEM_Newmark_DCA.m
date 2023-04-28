%----------------------------- DCA -------------------------------%
function [disp,vel,acc,M,C,K,fem] = FEM_Newmark_DCA(fem,opt,E,V,nb)
% Initial Effective Stiffness Matric
if ~isempty(fem.Kn0)
    Kn0 = fem.Kn0;
else
    fem.Kn0 = IniKn(fem,opt);
    Kn0 = fem.Kn0;
end
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
disp = zeros(2*fem.NNode,fem.NStep+1); vel = disp; acc = disp;
disp(:,1) = fem.u0; vel(:,1) = fem.v0; % Initial conditions
[M,C,K] = GlobalMCK(fem,E,V); % Compute stiffness and mass matrices
if ~isempty(fem.Mass)
  DOF = [2*fem.Mass(:,1)-1;2*fem.Mass(:,1)];
  M(DOF,DOF)= M(DOF,DOF) + [fem.Mass(:,2);fem.Mass(:,2)];
end
if ~isempty(fem.ag)
  fem.Fa = -M*ones(2*fem.NNode,1).*fem.ag;
  Ft = fem.Fext+fem.Fa;
else
  Ft = fem.Fext;
end
Kn = K+a0*M+a1*C;
% Threshold 1 of adaptively construct basis vectors
DeltK = (Kn-Kn0)./max(max(Kn0)); MaxDeltK = full(max(max(DeltK)));
% MaxKErr = 1;
acc(f,1) = M(f,f)\(Ft(f,1)-K(f,f)*disp(f,1)-C(f,f)*vel(f,1));
for It = 2:fem.NStep+1
    Fn = Ft(:,It)+M*(a0*disp(:,It-1)+a2*vel(:,It-1)+a3*acc(:,It-1))...
        +C*(a1*disp(:,It-1)+a4*vel(:,It-1)+a5*acc(:,It-1));
    if It == 2
        if MaxDeltK > 0.02
            deltK = Kn-Kn0;
            rb = zeros(2*fem.NNode,nb);
            [rb(f,1),L,s] = SolveLinSys(Kn0(f,f),Fn(f));
            for irb = 2:nb
                vec1 = deltK(f,f)*rb(f,irb-1);
                vec2(s) = L'\(L\vec1(s));
                rb(f,irb) = -vec2;
            end
            [rb,~,~] = svd(rb,0);
        else
            rb = fem.rb;
        end
        Kr = rb'*Kn*rb;
        % CA
        Fr = rb'*Fn;
        z = Kr\Fr;        
        disp(f,It) = rb(f,:)*z;
        Fn0 = Fn;
    else
        % Threshold 2 of adaptively construct basis vectors        
        DeltFn = Fn-Fn0; MaxDeltFn1 = full(max(abs(DeltFn)));
        MaxDeltFn1 = 1;
        if MaxDeltFn1 > 1
            deltK = Kn-Kn0;
            nb = 6;
            rb = zeros(2*fem.NNode,nb);
            [rb(f,1),L,s] = SolveLinSys(Kn0(f,f),Fn(f));
            for irb = 2:nb
                vec1 = deltK(f,f)*rb(f,irb-1);
                vec2(s) = L'\(L\vec1(s));
                rb(f,irb) = -vec2;
            end
            rb = orth(rb);
            Kr = rb'*Kn*rb;
            % CA
            Fr = rb'*Fn;
            z = Kr\Fr;
            disp(f,It) = rb(f,:)*z;
            Fn0 = Fn;
        else
            Fr = rb'*Fn;
            z = Kr\Fr;
            disp(f,It) = rb(f,:)*z;
            Fn0 = Fn;
        end
    end    
    acc(:,It) = a0*(disp(:,It)-disp(:,It-1))-a2*vel(:,It-1)-a3*acc(:,It-1);
    vel(:,It) = vel(:,It-1)+a6*acc(:,It-1)+a7*acc(:,It);
end
fem.Kn0 = Kn; fem.rb = rb;
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
%% ----------------------- INITIAL EFFECTIVE STIFFNESS MATRICES
function [Kn] = IniKn(fem,opt)
P=opt.P; z = ones(size(P,2),1);
[E,~,V,~] = opt.MatIntFnc(P*z);
[M,C,K] = GlobalMCK(fem,E,V);
% Parameter
gama = fem.gamma; beta = fem.beta; dt = fem.Tmax/fem.NStep+1; 
a0 = 1/beta/dt^2;
a1 = gama/beta/dt;
Kn = K+a0*M+a1*C;
%-------------------------------------------------------------------------%