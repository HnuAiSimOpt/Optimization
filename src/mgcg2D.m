% FUNCTION mgcg - MULTIGRID PRECONDITIONED CONJUGATE GRADIENTS
function [i,relres,u] = mgcg2D(A,b,u,Lfac,Ufac,Pu,nl,nswp,tol,maxiter)
r = b - A{1,1}*u;
res0 = norm(b);
% Jacobi smoother
omega = 0.8;
invD = cell(nl-1,1);
for l = 1:nl-1
    invD{l,1}= 1./spdiags(A{l,1},0);
end
for i = 1:1e6 
    z = VCycle(A,r,Lfac,Ufac,Pu,1,nl,invD,omega,nswp);
    rho = r'*z;
    if i == 1
        p = z;
    else
        beta = rho/rho_p;
        p = beta*p + z;
    end
    q = A{1,1}*p;
    dpr = p'*q;
    alpha = rho/dpr;
    u = u + alpha*p;
    r = r - alpha*q;
    rho_p = rho;
    relres = norm(r)/res0;
    if relres < tol || i>=maxiter
        break
    end
end
end


% FUNCTION VCycle - COARSE GRID CORRECTION
function z = VCycle(A,r,Lfac,Ufac,Pu,l,nl,invD,omega,nswp)
z = 0*r;
z = smthdmpjac(z,A{l,1},r,invD{l,1},omega,nswp);
Az = A{l,1}*z;
d = r - Az;
dh2 = Pu{l,1}'*d;
if (nl == l+1)
    vh2 = Ufac \ (Lfac \ dh2);
else
    vh2 = VCycle(A,dh2,Lfac,Ufac,Pu,l+1,nl,invD,omega,nswp);
end
v = Pu{l,1}*vh2;
z = z + v;
z = smthdmpjac(z,A{l,1},r,invD{l,1},omega,nswp);
end


% FUNCTIODN smthdmpjac - DAMPED JACOBI SMOOTHER
function [u] = smthdmpjac(u,A,b,invD,omega,nswp)
for i = 1:nswp
    u = u - omega*invD.*(A*u) + omega*invD.*b;
end
end


