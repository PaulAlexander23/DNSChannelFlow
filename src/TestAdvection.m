%--------------------------------------
% Driver routine for cavity flow
%--------------------------------------

%clear all
close all

global Rd Ad Pd
global Rp Ap Pp
global xLen yLen
global dx dy
global N M
global dt Re
global q2D

%...domain size
xLen  = 2;
yLen  = 1;

N     = 128;
M     = 64;

x       = linspace(0,xLen,N+1);
y       = linspace(0,yLen,M+1);
xc      = (x(1:end-1)+x(2:end))/2;
yc      = (y(1:end-1)+y(2:end))/2;
[yy,xx] = meshgrid(yc,xc);
dx      = xLen/N;
dy      = yLen/M;
dt      = min(dx,dy)/1.5;

q2   = zeros(N,M);
q2(96:128,48:52) = 1;
ax    = ones(N,M)*1;
ay    = ones(N,M)*0.1;

q2 = q2(:);
ax  = ax(:);
ay  = ay(:);

q2S =  zeros(N,1);
q2N =  zeros(N,1);
q2W =  zeros(1,M);
q2E =  zeros(1,M);

for i=1:50
    
    %...semi-Lagrangian advection
    q2    = SemiLagrAdvect(ax,ay,q2,q2S,q2N,q2W,q2E);
    
    if (mod(i,1)==0)
        k = 2;
        fprintf('time step %i \n',i)
        %...graphics output
        q2D = reshape(q2,N,M);
        figure(1);
        surf(xx,yy,q2D,'EdgeColor','None')
        axis image
        axis([0 xLen 0 yLen]);
        drawnow
    end
end