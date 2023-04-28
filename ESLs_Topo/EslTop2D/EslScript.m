%------------------------------ EslScript 2D-----------------------------%
% ''An efficient online successive reanalysis                            %
%   method for dynamic topology optimization''                           %
%------------------------------------------------------------------------%
clear; clc; close all
restoredefaultpath; addpath(genpath('./'));% Use all folders and subfolders
set(0,'defaulttextinterpreter','latex')
% Finite Element Model
[Node,Element,Supp,Load] = PolyMesher(@EyeBarDomain,10000,20);
Tmax = 0.05; % Simulation time
NStep = (size(Load,2)-3)/2; % Number of time steps
O = zeros(2*size(Node,1),1);
%% Initialization FEM
fem = struct(...
  'NNode',size(Node,1),...     % Number of nodes
  'NElem',size(Element,1),...  % Number of elements
  'Node',Node,...              % [NNode x 3] array of nodes
  'Element',{Element},...      % [NElement x Var] cell array of elements
  'Supp',Supp,...              % Array of supports
  'Load',Load,...              % Array of loads
  'Mass',[],...                % Array of lumped masses
  'Ar',[10,1e-5],...           % Rayleigh damping param. C=Ar(1)*M+Ar(2)*K
  'u0',O,...                   % Initial displacement vector
  'v0',O,...                   % Initial velocity vector
  'Thickness',0.01,...         % Element thickness
  'Tmax',Tmax,...              % Simulation time
  'NStep',NStep, ...           % Total number of steps 
  'beta', 0.25, ...            % beta parameter for Newmark-beta method
  'gamma',0.5, ...             % gamma parameter for Newmark-beta method
  'NIter',0,...                % Total number of optimization iterations  
  'LL',[],...                  % Vector of DOF index for U_DOF objective
  'Nu0',0.3,...                % Poisson's ratio of solid material
  'E0',2.1e5,...               % Young's modulus of solid material
  'rho',7.9e-9,...             % Mass density of solid material
  'ag',[],...                  % Ground acceleration
  'Kn0',[],...                 % Initial effective stiffness matric
  'rb',[],...                  % Basis vectors  
  'Reg',0 ...                  % Tag for regular meshes  
   );
clear Node Element Supp Load 
%% Initialization OPT
R = 0.1; q = 1; VolFrac = 0.5;
m = @(y)MatIntFnc(y,'SIMP',3);
fem.SElem = []; % Elements in passive solid regions
ElemInd{1} = (1:fem.NElem)'; % Element indices for volume constraint j
P = PolyFilter(fem,R,q,'X');% Filter
zIni = ones(size(P,2),1); 
for ii=1:length(VolFrac); zIni(ElemInd{ii})=VolFrac(ii); end % Initial DVs
opt = struct(...               
  'zMin',0.0,...               % Lower bound for design variables
  'zMax',1.0,...               % Upper bound for design variables
  'zIni',zIni,...              % Initial design variables
  'MatIntFnc',m,...            % Handle to material interpolation fnc.
  'P',P,...                    % Matrix that maps design to element vars.
  'VolFrac',VolFrac,...        % Specified volume fraction cosntraint
  'Tol',0.01,...               % Convergence tolerance on design vars.  
  'MaxIter',20,...             % Max. number of optimization iterations
  'OCMove',0.2,...             % Allowable move step in OC update scheme
  'OCEta',0.5 ...              % Exponent used in OC update scheme
   );
clear R q VolFrac m P zIni
%% Optimization
fem = PreComputations(fem);
figure; t = cputime;
B = 1; p_i = 0:1.5:9; 
for ii=1:length(p_i) % Continuation on the RAMP penalty parameter
  if ii==length(p_i); opt.MaxIter = 100; end
  disp(['current p: ', num2str(p_i(ii)), '  current B: ', num2str(B)]);
  opt.MatIntFnc = @(y)MatIntFnc(y,'RAMP-H1',[p_i(ii),B,0.5]);
  [opt.zIni,V,fem] = EslTop(fem,opt);
  B = min(B+2,10);
end
t0 = cputime-t;
% Total optimization time
if t0<=60, fprintf('Optimization time: %i seconds \n', round(t0))
elseif t0<=3600, fprintf('Optimization time: %1.1f minutes \n', t0/60) 
elseif t0<=86400, fprintf('Optimization time: %1.1f hours \n', t0/3600) 
else, fprintf('Optimization time: %1.1f days \n', t0/86400) 
end