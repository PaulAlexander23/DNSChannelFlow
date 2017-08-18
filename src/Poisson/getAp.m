function getAp(N,M)
%======================================
% set up a hierarchy of Laplace
% operators (as cell arrays)
%======================================

global Rp Ap Pp
global dx dy

%...number of levels and initialization
kk   = length(Pp);
Ap   = cell(kk+1,2);

%...Laplace operator
AAx      = spdiags(ones(N,1)*[1 -2 1]/dx/dx,-1:1,N,N);
AAx(1,1) = -1/dx/dx;
AAx(N,N) = -3/dx/dx;

AAy      = spdiags(ones(M,1)*[1 -2 1]/dy/dy,-1:1,M,M);
AAy(1,1) = -1/dy/dy;
AAy(M,M) = -1/dy/dy;

AAA  = kron(speye(M),AAx) +kron(AAy,speye(N));
AAA(3*N/4+1:N,3*N/4+1:N) = AAA(3*N/4+1:N,3*N/4+1:N) - 2/dy/dy*speye(N/4);
AAA(1,1) = -3/dx/dx;
%...lower-level operators (Galerkin condition)
Ap{1,1} = AAA;
for i=2:kk+1
    Ap{i,1} = Rp{i-1,1}*Ap{i-1,1}*Pp{i-1,1};
end

%...Laplace operator
AAx(N,N) = -1/dx/dx;

AAA  = kron(speye(M),AAx) + kron(AAy,speye(N));
AAA(N/2+1:N,N/2+1:N) = AAA(N/2+1:N,N/2+1:N) - 2/dy/dy*speye(N/2);
%...lower-level operators (Galerkin condition)
Ap{1,2} = AAA;
for i=2:kk+1
    Ap{i,2} = Rp{i-1,2}*Ap{i-1,2}*Pp{i-1,2};
end
end