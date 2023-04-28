function [u] = mgcg(A,b,u,Pu,nl,nswp,tol,maxiter)
r = b - A{1,1}*u;
res0 = norm(b);
% Jacobi smoother
omega = 0.6;
invD = cell(nl-1,1);
for l = 1:nl-1
    invD{l,1}= 1./spdiags(A{l,1},0);
end
for i = 1:1e6 
    z = VCycle(A,r,Pu,1,nl,invD,omega,nswp);
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
%fprintf('cgiters:%i,  relres:%f\n',i,relres);
end