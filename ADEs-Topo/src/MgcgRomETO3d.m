function [TFE,Comp,err,loop,Tiji] = MgcgRomETO3d(nelx,nely,nelz,volfrac,er,rmin,start)
% nelx = 60; nely = 20; nelz = 12; volfrac = 0.3; er = 0.04; rmin = 1.5; start = 25;
% 初始材料参数设置
E0 = 1;   dotmin = 1e-9;   Emin = E0*dotmin;   nu = 0.3;   penal = 3.0;
vx = ones(nely,nelx,nelz);
% 初始算法参数设置
loopMax = 600;   stop = 1e-6;   showflag = 0;
Nb = 4;  nl=3;  epsilon1 = 0.04;  epsilon2 = 0.06;
nswp = 1;  cgtol = 1e-6;  cgmax = 50;   ngrid = 10;

% 延拓操作
Pu = cell(nl-1,1);
for l = 1:nl-1
    [Pu{l,1}] = prepcoarse(nelx/2^(l-1),nely/2^(l-1),nelz/2^(l-1));
end

% 载荷自由度定义
% il = nelx/2; jl = nely/2; kl = nelz;                    % Coordinates
il = ones(1,nely+1)*(nelx/2);
jl = 0:nely;
kl = ones(1,nely+1)*nelz;
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:);                                 % DOFs 

% 固支自由度定义
if1 = [0;0;nelx;nelx];
jf1 = [nely;0;nely;0];
kf1 = [0;0;0;0];
fixednid = kf1*(nelx+1)*(nely+1)+if1*(nely+1)+(nely+1-jf1);    % Node IDs
fixeddofs = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs  


% FEA准备
KE = KE3d(nu);
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-1,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddofs);

% Null space elimination of supports
N = ones(ndof,1); N(fixeddofs) = 0; Null = spdiags(N,0,ndof,ndof);

% 单元节点和编号
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);

% 节点坐标
nodeidspro = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidzpro = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeidspro = repmat(nodeidspro,size(nodeidzpro))+repmat(nodeidzpro,size(nodeidspro));
elenod = repmat(nodeidspro(:),1,8)+repmat([1 nely+[2 1] 0 (nely+1)*(nelx+1)+[1 nely+[2 1] 0]],nelx*nely*nelz,1);
[nodex, nodey, nodez] = meshgrid(0:1:nelx,nely:-1:0,0:1:nelz);

%过滤准备
[H,Hs] = TopFilter(nelx,nely,nelz,rmin);

% 投影准备
[s,t,w] = meshgrid(-1+1/ngrid:2/ngrid:1-1/ngrid,-1+1/ngrid:2/ngrid:1-1/ngrid,-1+1/ngrid:2/ngrid:1-1/ngrid);
N1 = (1 - s(:)).*(1 - t(:)).*(1 - w(:))/8;
N2 = (1 + s(:)).*(1 - t(:)).*(1 - w(:))/8;
N3 = (1 + s(:)).*(1 + t(:)).*(1 - w(:))/8;
N4 = (1 - s(:)).*(1 + t(:)).*(1 - w(:))/8;
N5 = (1 - s(:)).*(1 - t(:)).*(1 + w(:))/8;
N6 = (1 + s(:)).*(1 - t(:)).*(1 + w(:))/8;
N7 = (1 + s(:)).*(1 + t(:)).*(1 + w(:))/8;
N8 = (1 - s(:)).*(1 + t(:)).*(1 + w(:))/8;

% 迭代初始化
c = zeros(loopMax,1);
vol = 1.0; loop = 0; change = 1.0;
Erb = zeros(1,loopMax);
fullnum = 0;  e1num = 0;  FT = 0;
Sizefree = size(freedofs,2);  RB = zeros(Sizefree,Nb);
% 开始迭代
T = zeros(1,loopMax);
Tiji = zeros(1,loopMax);
while change > stop && loop < loopMax
  loop = loop + 1;  vol = max(vol*(1-er),volfrac);
  if loop > 2
    olddcnd_2 = olddcnd_1; olddcnd_1 = dcnd;
  end
  sK = reshape(KE(:)*(vx(:)'*E0+(1-vx(:))'*Emin),576*nele,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  tic;
  if loop < start
    Kmg = cell(nl,1);
    Kmg{1,1} = Null'*K*Null - (Null-speye(ndof,ndof));
    for l = 1:nl-1
      Kmg{l+1,1} = Pu{l,1}'*(Kmg{l,1}*Pu{l,1});
    end
    U = mgcg(Kmg,F,U,Pu,nl,nswp,cgtol,cgmax);
    U(fixeddofs) = 0;
    fullnum = fullnum + 1;    Erb(loop) = nan;
  elseif loop == start
    Kmg = cell(nl,1);
    Kmg{1,1} = Null'*K*Null - (Null-speye(ndof,ndof));
    for l = 1:nl-1
      Kmg{l+1,1} = Pu{l,1}'*(Kmg{l,1}*Pu{l,1});
    end
    U = mgcg(Kmg,F,U,Pu,nl,nswp,cgtol,cgmax);
    U(fixeddofs) = 0;
    norm2 = norm(U(freedofs));   RB(:,1) = U(freedofs)/norm2;
    fullnum = fullnum + 1;    Erb(loop) = nan; 
  elseif (loop > start) && (loop <= start+Nb)
    Kmg = cell(nl,1);
    Kmg{1,1} = Null'*K*Null - (Null-speye(ndof,ndof));
    for l = 1:nl-1
      Kmg{l+1,1} = Pu{l,1}'*(Kmg{l,1}*Pu{l,1});
    end
    U = mgcg(Kmg,F,U,Pu,nl,nswp,cgtol,cgmax);
    U(fixeddofs) = 0;
    U0 = U(freedofs);
    for i = 1:loop-start
      aij = U(freedofs)'*RB(:,i);
      U0 = U0 - aij*RB(:,i);
    end
    norm2 = norm(U0);  RB(:,loop-start+1) = U0/norm2;
    fullnum = fullnum + 1;    Erb(loop) = nan;
  elseif loop > start+Nb
    switch FT
      case 0
        Kr = RB'*K(freedofs,freedofs)*RB;   Fr = RB'*F(freedofs);
        Uz = Kr\Fr;   U0 = RB*Uz;
        norm_err = norm(K(freedofs,freedofs)*U0-F(freedofs));
        norm_F = norm(F(freedofs));
        Erb(loop) = norm_err/norm_F;
        e1num = e1num +1;
        if Erb(loop) < epsilon1
          U(freedofs) = U0; FT = 0;
        elseif (Erb(loop) >= epsilon1) && (Erb(loop) < epsilon2)
          U(freedofs) = U0; FT = 1;
        else
          RB(:,1:end-1) = RB(:,2:end);
          Kmg = cell(nl,1);
          Kmg{1,1} = Null'*K*Null - (Null-speye(ndof,ndof));
          for l = 1:nl-1
            Kmg{l+1,1} = Pu{l,1}'*(Kmg{l,1}*Pu{l,1});
          end
          U = mgcg(Kmg,F,U,Pu,nl,nswp,cgtol,cgmax);
          U(fixeddofs) = 0;
          U0 = U(freedofs);
          fullnum = fullnum + 1;
          for i = 1:Nb-1
            aij = U(freedofs)'*RB(:,i);
            U0 = U0 - aij*RB(:,i);
          end
          norm2 = norm(U0);   RB(:,end) = U0/norm2;
          FT = 0;
        end
      case 1
        RB(:,1:end-1) = RB(:,2:end);
        Kmg = cell(nl,1);
        Kmg{1,1} = Null'*K*Null - (Null-speye(ndof,ndof));
        for l = 1:nl-1
          Kmg{l+1,1} = Pu{l,1}'*(Kmg{l,1}*Pu{l,1});
        end
        U = mgcg(Kmg,F,U,Pu,nl,nswp,cgtol,cgmax);
        U(fixeddofs) = 0;
        U0 = U(freedofs);
        fullnum = fullnum + 1;
        for i = 1:Nb-1
          aij = U(freedofs)'*RB(:,i);
          U0 = U0 - aij*RB(:,i);
        end
        norm2 = norm(U0);   RB(:,end) = U0/norm2;
        FT = 0;
    end
  end
  T(loop) = toc;
  % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
  c(loop) = sum(sum(sum((vx.*E0+(1-vx).*Emin).*ce)));
  dc = ((1-vx)*dotmin.^(penal-1)+vx)*E0.*ce;
  % 过滤
  dcnd = reshape(H*dc(:)./Hs,nely+1,nelx+1,nelz+1);
  % 保证进化过程稳定
  if loop == 1
    olddcnd_2 = dcnd;
  elseif loop == 2
    olddcnd_1 = dcnd;
  else
    dcnd = (dcnd+olddcnd_1+olddcnd_2)/3.;
  end
  % ETO 设计变量更新
  L1 = min(dcnd(:)); L2 = max(dcnd(:));
  Ls = (L1+L2)/2.0;
  while (L2-L1)/abs(L1+L2) > 1.0e-9
    Ls = (L1+L2)/2.0;
    dcth = dcnd(:)-Ls;
    MaxDcth = max(dcth(elenod(:,:)),[],2);
    MinDcth = min(dcth(elenod(:,:)),[],2);
    index1 = find(MaxDcth<0);   index2 = find(MinDcth>0);
    vx(index1) = 0;  vx(index2) = 1;
    index3 = setdiff([1:1:nelz*nely*nelx],union(index1,index2));
    [~,index3Size] = size(index3);
    for j = 1:index3Size
      i = index3(j);
      ps = N1 * dcth(elenod(i,1)) + N2*dcth(elenod(i,2)) + N3*dcth(elenod(i,3)) + N4*dcth(elenod(i,4))+...
        N5 * dcth(elenod(i,5)) + N6*dcth(elenod(i,6)) + N7*dcth(elenod(i,7)) + N8*dcth(elenod(i,8));
      vx(i) = length(find( ps >= 0 ))/length(s(:));
    end
    if mean(vx(:)) - vol > 0
      L1 = Ls;
    else
      L2 = Ls;
    end
  end
  % PRINT RESULTS
  Tiji(loop) = sum(sum(sum(vx)))/(nelx*nely*nelz);
  if loop > 10; change =abs(sum(c(loop-9:loop-5))-sum(c(loop-4:loop)))/sum(c(loop-4:loop)); end
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c(loop)) ...
    ' Vol.: ' sprintf('%6.3f',sum(sum(sum(vx)))/(nelx*nely*nelz)) ...
    ' ch.: ' sprintf('%6.3e',change )])
  % PLOT DENSITIES
  if showflag
    figure(1)
    clf;
    p = patch(isocaps(nodex, nodey, nodez,dcnd,Ls),'FaceColor','blue','EdgeColor','none');
    q = patch(isosurface(nodex, nodey, nodez,dcnd,Ls),'FaceColor','red','EdgeColor','none');
    isonormals(nodex, nodey, nodez,dcnd,p);
    isonormals(nodex, nodey, nodez,dcnd,q);
    %rotate(p,[1 0 0],90);  rotate(q,[1 0 0],90);
    daspect([1 1 1]); 
    view(37.5,30); axis tight; axis off; camlight;
  end
end
TFE = sum(T(:));
Comp = c(1:loop,1);
err = Erb;
figure(1)
clf;
p = patch(isocaps(nodex, nodey, nodez,dcnd,Ls),'FaceColor','blue','EdgeColor','none');
q = patch(isosurface(nodex, nodey, nodez,dcnd,Ls),'FaceColor','red','EdgeColor','none');
isonormals(nodex, nodey, nodez,dcnd,p);
isonormals(nodex, nodey, nodez,dcnd,q);
%rotate(p,[1 0 0],90);  rotate(q,[1 0 0],90);
daspect([1 1 1]);
view(45,15); axis tight; axis off; camlight;
figure(3)
clf;
p = patch(isocaps(nodex, nodey, nodez,dcnd,Ls),'FaceColor','blue','EdgeColor','none');
q = patch(isosurface(nodex, nodey, nodez,dcnd,Ls),'FaceColor','red','EdgeColor','none');
isonormals(nodex, nodey, nodez,dcnd,p);
isonormals(nodex, nodey, nodez,dcnd,q);
%rotate(p,[1 0 0],90);  rotate(q,[1 0 0],90);
daspect([1 1 1]);
view(80,15); axis tight; axis off; camlight;
