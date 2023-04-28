function [T,Comp,Res,mgcgVec] = comp2DbyCA(nelx,nely,vol,rmin,nROM,nCA,tol,nl,On)
% 默认输入参数
cgtol = 1e-6;  cgmax = 200;
% 材料参数及算法参数
E0 = 1;  E_factor = 1E-9;  Emin = E0*E_factor;  nu = 0.3;   penal = 3;
Maxiter = 150;  Exit_tol = 1e-20;  iters = 0;  change = 1.;  start = On-nROM+1;

% Prologation operators
Pu = cell(nl-1,1); 
for l = 1:nl-1  
    [Pu{l,1}] = Prepcoarse2D(nely/2^(l-1),nelx/2^(l-1)); 
end

% 有限元分析准备
[KE,iK,jK,edofMat] = Perp2D(nu,nelx,nely);

% 定义边界条件
F = sparse(2*(nely+1)*(nelx/2+1),1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union(2*nely+1:2*(nely+1), 2*(nely+1)*(nelx+1));
alldofs = (1:2*(nely+1)*(nelx+1));
freedofs = setdiff(alldofs,fixeddofs);
nfree = size(freedofs,2);

% Null space elimination of supports
ndof = 2*(nelx+1)*(nely+1);
N = ones(ndof,1); N(fixeddofs) = 0; Null = spdiags(N,0,ndof,ndof);

% 过滤准备
[H,Hs] = Filter2D(nelx,nely,rmin);   % 过滤准备

% 初始化迭代
x = repmat(vol,nely,nelx);   xPhys = x;   ROM = zeros(nfree,nROM);
norm_F = norm(F(freedofs));
Res = zeros(Maxiter,1);   Comp = zeros(Maxiter,1); T = zeros(Maxiter,1);
mgcgVec = zeros(Maxiter,1);
while change>Exit_tol && iters<Maxiter
    iters = iters + 1;
    fprintf(' It.:%5i ',iters);
    if iters < start
        sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        Kk = sparse(iK,jK,sK);  Kk = (Kk+Kk')/2;
        
        tic;
        K = cell(nl,1); K{1,1} = Null'*Kk*Null - (Null-speye(ndof,ndof));
        for l = 1:nl-1; K{l+1,1} = Pu{l,1}'*(K{l,1}*Pu{l,1}); end
        Lfac = chol(K{nl,1},'lower'); Ufac = Lfac';
        [mgcgVec(iters),~,U] = mgcg2D(K,F,U,Lfac,Ufac,Pu,nl,1,cgtol,cgmax);  U(fixeddofs)=0;
        Res(iters) = nan;
        T(iters) = toc;
    else
        oldK = Kk;
        sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        Kk = sparse(iK,jK,sK);  Kk = (Kk+Kk')/2;

        tic;
        if iters == start
            K = cell(nl,1);  K{1,1} = Null'*Kk*Null - (Null-speye(ndof,ndof));
            for l = 1:nl-1; K{l+1,1} = Pu{l,1}'*(K{l,1}*Pu{l,1}); end
            Lfac = chol(K{nl,1},'lower');  Ufac = Lfac';
            [mgcgVec(iters),~,U] = mgcg2D(K,F,U,Lfac,Ufac,Pu,nl,1,cgtol,cgmax);  U(fixeddofs)=0;
            ROM(:,iters-start+1) = U(freedofs)/norm(U(freedofs));
            Res(iters) = nan;
        elseif (iters>start) && (iters-start<nROM)
            K = cell(nl,1);  K{1,1} = Null'*Kk*Null - (Null-speye(ndof,ndof));
            for l = 1:nl-1; K{l+1,1} = Pu{l,1}'*(K{l,1}*Pu{l,1}); end
            Lfac = chol(K{nl,1},'lower');  Ufac = Lfac';
            [mgcgVec(iters),~,U] = mgcg2D(K,F,U,Lfac,Ufac,Pu,nl,1,cgtol,cgmax);  U(fixeddofs)=0;
            ROM(:,iters-start+1) = U(freedofs)/norm(U(freedofs));
            Res(iters) = nan;
        else
            cab = CA2D(ROM,Kk(freedofs,freedofs),oldK(freedofs,freedofs),U(freedofs),nfree,nCA);
            Kr = cab'*Kk(freedofs,freedofs)*cab;   Fr = cab'*F(freedofs);
            U0 = cab*(Kr\Fr);
            norm_e = norm(Kk(freedofs,freedofs)*U0-F(freedofs));
            Res(iters) = norm_e/norm_F;
            if Res(iters) < tol
                U(freedofs) = U0;
            else
                K = cell(nl,1);  K{1,1} = Null'*Kk*Null - (Null-speye(ndof,ndof));
                for l = 1:nl-1; K{l+1,1} = Pu{l,1}'*(K{l,1}*Pu{l,1}); end
                Lfac = chol(K{nl,1},'lower');  Ufac = Lfac';
                [mgcgVec(iters),~,U] = mgcg2D(K,F,U,Lfac,Ufac,Pu,nl,1,cgtol,cgmax);  U(fixeddofs)=0;
                ROM(:,1:end-1) = ROM(:,2:end);
                ROM(:,end) = U(freedofs)/norm(U(freedofs));
            end
        end
        T(iters) = toc;
    end
    % 目标函数和敏度
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    Comp(iters) = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx);
    
    % 过滤/敏度修正
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
    
    % OC优化
    l1 = 0; l2 = 1e9;  move = 0.1;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
        xPhys(:) = (H*xnew(:))./Hs;
        if sum(xPhys(:)) > vol*nelx*nely, l1 = lmid; else l2 = lmid; end
    end
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
  % PRINT RESULTS
  fprintf(' Obj.:%11.4f Vol.:%7.3f ch.:%7.3f err.:%5f rmin.:%.2f\n', ...
      Comp(iters),mean(xPhys(:)),change, Res(iters),rmin);
end
% PLOT DENSITIES
figure(1)
colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
figure(2)
semilogy([1:iters],Res(1:iters),[1:iters],ones(1,iters)*tol,'r');
