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
    figure(2)
    subplot(1,2,1); [Fig] = InitialPlot(fem,opt,opt.PF*z,[30,30]);
    subplot(1,2,2); [Fig1] = InitialPlot(fem,opt,opt.PF*z,[0,90]);
    Iter = 0;
    MaxIterS = 30; % Max. number of static optimization iterations
    while (Iter<MaxIterS) && (Change>Tol)
        tic
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
        % Outputda results
        NIter = fem.NIter;
        fprintf('It: %i \t Objective: %1.3f\tChange: %1.3f\n',NIter,f,Change);
        subplot(1,2,1);
        ZPlot=reshape(Fig.PPlot*opt.PF*z,size(Fig.X));
        [faces,verts] = isosurface(Fig.X,Fig.Y,Fig.Z,ZPlot,0.5);
        set(Fig.Handle,'Faces', faces, 'Vertices', verts);
        subplot(1,2,2);
        ZPlot=reshape(Fig1.PPlot*opt.PF*z,size(Fig1.X));
        [faces,verts] = isosurface(Fig1.X,Fig1.Y,Fig1.Z,ZPlot,0.5);
        set(Fig1.Handle,'Faces', faces, 'Vertices', verts);
        drawnow
    end
end
%------------------------------------------------------- OBJECTIVE FUNCTION
function [f,dfdE,dfdV,vem] = ObjectiveFnc(vem,E,F,V)
K = sparse(vem.i,vem.j,E(vem.e).*vem.k0);
K = (K+K')/2;
n = size(F,2);
U = zeros(3*vem.NNode,n);
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
function [g,dgdE,dgdV,vem] = ConstraintFnc(vem,E,V,VolFrac)
g = sum(vem.ElemVolume.*V)/sum(vem.ElemVolume)-VolFrac;
dgdE = zeros(size(E));
dgdV = vem.ElemVolume/sum(vem.ElemVolume);
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
function [Fig] = InitialPlot(vem,opt,z0,shijiao)
BdBox = vem.Domain('BdBox');
dmin = 0.5*(sum(vem.ElemVolume)/vem.NElem)^(1/3);
[Fig.X,Fig.Y,Fig.Z]=...
  meshgrid(BdBox(1)-5*dmin/2:dmin:BdBox(2)+5*dmin/2,...
           BdBox(3)-5*dmin/2:dmin:BdBox(4)+5*dmin/2,...
           BdBox(5)-5*dmin/2:dmin:BdBox(6)+5*dmin/2);
h_E=(sum(vem.ElemVolume)/length(vem.Element))^(1/3);
PS1=[Fig.X(:),Fig.Y(:),Fig.Z(:)];
% Dist = vem.Domain('Dist',PS1); I = find(Dist(:,end)<0); PS1=PS1(I,:);
PS2=[opt.NodeD(:,1),opt.NodeD(:,2),opt.NodeD(:,3)];
Fig.PPlot = PolyFilter3D(PS1,PS2,h_E,2);
ZPlot=reshape(Fig.PPlot*opt.PF*z0,size(Fig.X));
[faces,verts] = isosurface(Fig.X,Fig.Y,Fig.Z,ZPlot,0.5);
cla, hold on, view(shijiao), rotate3d on, axis equal
box
axis off
Fig.Handle=patch('Faces', faces, 'Vertices', verts,...
     'FaceColor', 'r', 'EdgeColor', 'none',...
     'FaceLighting','gouraud','AmbientStrength',0.5);
hold off;
camlight
drawnow