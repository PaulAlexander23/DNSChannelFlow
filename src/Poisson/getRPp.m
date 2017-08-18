function getRPp(N,M)
%======================================
% set up a hierarchy of restriction
% and prolongation matrices (as cell
% arrays) (Neumann version)
%======================================

global Rp Ap Pp
global ilevmin

%...number of levels and initialization
kk      = log(N)/log(2)-1;
ll      = log(M)/log(2)-1;
kl      = min(kk,ll);
ilevmin = kl;
Rp      = cell(kl,2);
Pp      = cell(kl,2);

%...set up prolongation
NN = N/2;
MM = M/2;
for i=1:kl
    
    %Neumann East, Dirichlet West
    PPx = sparse(2*NN,NN);
    for j=1:NN-1
        PPx(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
    end
    PPx(1,1)     = 1;
    PPx(2*NN,NN) = 0.5;
    
    %Neumann North and South
    PPy1 = sparse(2*MM,MM);
    for j=1:MM-1
        PPy1(2*j:2*j+1,j:j+1) = [0.75 0.25; 0.25 0.75];
    end
    PPy1(1,1)     = 1;
    PPy1(2*MM,MM) = 1;
    
    %Neumann North, Dirichlet South
    PPy2 = PPy1;
    PPy2(1,1)     = 0.5;
    
    %Splicing for Omega 1
    E = cat(2,eye(3*NN/4), zeros(3*NN/4,NN/4));
    F = kron(eye(MM),E);
    G = zeros(NN/4,NN);
    G(1:NN/4,3*NN/4+1:end) = eye(NN/4);
    H = kron(eye(MM),G);
    
    Pp{i,1} = kron(PPy1,PPx(:,1:3*NN/4))*F + kron(PPy2,PPx(:,3*NN/4+1:NN))*H;
    
    %Adding a single Dirichlet point
    Pp{i,1}(1,1) = 0.5;
    
    %Changing East boundary to Neumann
    PPx(2*NN,NN) = 1;
    
    %Splicing for Omega 2
    E = cat(2,eye(NN/2), zeros(NN/2,NN/2));
    F = kron(eye(MM),E);
    G = zeros(NN/2,NN);
    G(1:NN/2,NN/2+1:end) = eye(NN/2);
    H = kron(eye(MM),G);
    
    Pp{i,2} = kron(PPy1,PPx(:,1:NN/2))*F + kron(PPy2,PPx(:,NN/2+1:NN))*H;
    
    NN = NN/2;
    MM = MM/2;
    
end

%...set up restriction (transpose of prolongation)
for i=1:kl
    Rp{i,1} = transpose(Pp{i,1});
    Rp{i,2} = transpose(Pp{i,2});
end
end

