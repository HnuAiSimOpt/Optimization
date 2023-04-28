function [z,V,fem] = EslTop(fem,opt)
Tol = opt.Tol*(opt.zMax-opt.zMin);
Change = 2*Tol; z = opt.zIni; P=opt.P;
[E,dEdy,V,dVdy] = opt.MatIntFnc(P*z);
IterDy = 0;
while (IterDy<opt.MaxIter) && (Change>Tol)
    IterDy = IterDy+1;
    % Run dynamics code
    [U,~,~,~,~,K,fem] = FEM_Newmark(fem,E,V);
    % Dynamic combined approximation (DCA)
    nb = 6; % Number of basis vectors 
    [U,~,~,~,~,K,fem] = FEM_Newmark_DCA(fem,opt,E,V,nb);     
    % Mode superposition method (MSM)
    nm = 1; % Number of modes
    [U,~,~,~,~,K,fem] = FEM_MSM(fem,E,V,nm);
    % Exact ESLs
    Feq = K*U;
    % Approximate ESLs by POD
    [U1,S,~] = svd(Feq,0); % SVD
    S = diag(S); % Singular value
    i = 0; Ratio = 0;
    while Ratio<0.9 % Given energy ratio
        i = i+1;
        Ratio = sum(S(1:i))/sum(S);
    end
    nPOD =  i; % Number of approximated ESLs
    Feq = U1(:,1:nPOD);
    figure(2); [FigHandle,FigData] = InitialPlot(fem,V);
    Iter = 0; 
    MaxIterS = 30; % Max. number of static optimization iterations     
    while (Iter<MaxIterS) && (Change>Tol)
        Iter = Iter + 1;
        fem.NIter = fem.NIter+1;
        % Compute cost functionals and analysis sensitivities
        [f,dfdE,dfdV,fem] = ObjectiveFnc(fem,E,Feq,V);
        [g,dgdE,dgdV,fem] = ConstraintFnc(fem,E,V,opt.VolFrac); 
        % Compute design sensitivities
        dfdz = P'*(dEdy.*dfdE + dVdy.*dfdV);
        dgdz = P'*(dEdy.*dgdE + dVdy.*dgdV);
        % Update design variable and analysis parameters
        [z,Change] = UpdateScheme(dfdz,g,dgdz,z,opt);
        [E,dEdy,V,dVdy] = opt.MatIntFnc(P*z);
        % Output results
        NIter = fem.NIter;
        fprintf('It: %i \t Objective: %1.3f\tChange: %1.3f\n',NIter,f,Change);
        figure(2); 
        set(FigHandle,'FaceColor','flat','CData',1-V(FigData));
        drawnow
    end
end
%------------------------------------------------------- OBJECTIVE FUNCTION
function [f,dfdE,dfdV,vem] = ObjectiveFnc(vem,E,F,V)
K = sparse(vem.i,vem.j,E(vem.e).*vem.k0);
K = (K+K')/2;
n = size(F,2);
U = zeros(2*vem.NNode,n);
U(vem.FreeDofs,:) = K(vem.FreeDofs,vem.FreeDofs)\F(vem.FreeDofs,:);
temp = 0;
for i = 1:n
    Ui = U(:,i);
    temp = temp + cumsum(-Ui(vem.i).*vem.k0.*Ui(vem.j));
end
temp = temp(cumsum(vem.ElemNDof.^2));
dfdE = [temp(1);temp(2:end)-temp(1:end-1)];
f = trace(F'*U); 
dfdV = zeros(size(V));
%------------------------------------------------------ CONSTRAINT FUNCTION
function [g,dgdE,dgdV,fem] = ConstraintFnc(fem,E,V,VolFrac)
if ~isfield(fem,'ElemArea')
  fem.ElemArea = zeros(fem.NElem,1);
  for el=1:fem.NElem
    vx=fem.Node(fem.Element{el},1); vy=fem.Node(fem.Element{el},2);
    fem.ElemArea(el) = 0.5*sum(vx.*vy([2:end 1])-vy.*vx([2:end 1]));
  end
end
g = sum(fem.ElemArea.*V)/sum(fem.ElemArea)-VolFrac;
dgdE = zeros(size(E));
dgdV = fem.ElemArea/sum(fem.ElemArea);
%----------------------------------------------- OPTIMALITY CRITERIA UPDATE
function [zNew,Change] = UpdateScheme(dfdz,g,dgdz,z0,opt)  
zMin=opt.zMin; zMax=opt.zMax;  
move=opt.OCMove*(zMax-zMin); eta=opt.OCEta;
l1=0; l2=1e6;  
while l2-l1 > 1e-4
  lmid = 0.5*(l1+l2);
  B = -(dfdz./dgdz)/lmid;
  zCnd = zMin+(z0-zMin).*B.^eta;
  zNew = max(max(min(min(zCnd,z0+move),zMax),z0-move),zMin);
  if (g+dgdz'*(zNew-z0)>0),  l1=lmid;
  else                       l2=lmid;  end
end
Change = max(abs(zNew-z0))/(zMax-zMin);
%------------------------------------------------------------- INITIAL PLOT
function [handle,map] = InitialPlot(fem,z0)
Tri = zeros(length([fem.Element{:}])-2*fem.NElem,3);
map = zeros(size(Tri,1),1); index=0;
for el = 1:fem.NElem
  for enode = 1:length(fem.Element{el})-2
    map(index+1) = el;
    Tri(index+1,:) = fem.Element{el}([1,enode+1,enode+2]);
    index = index + 1;
  end
end
handle = patch('Faces',Tri,'Vertices',fem.Node,'FaceVertexCData',...
               1-z0(map),'FaceColor','flat','EdgeColor','none');
axis equal; axis off; axis tight; colormap(gray); caxis([0 1]);