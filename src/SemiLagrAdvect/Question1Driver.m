%--------------------------------------
% Driver routine for Question one
%--------------------------------------

clear all
close all

global Rd Ad Pd
global Rp Ap Pp
global xLen yLen
global dx dy
global N M
global dt Re

%...Reynolds number
Re = 500;

%...domain size
xLen  = 0.8;
yLen  = 0.64;

N     = 80;
M     = 64;

x       = linspace(0,xLen,N+1);
y       = linspace(0,yLen,M+1);
xc      = (x(1:end-1)+x(2:end))/2;
yc      = (y(1:end-1)+y(2:end))/2;
[yy,xx] = meshgrid(yc,xc);
dx      = xLen/N;
dy      = yLen/M;
dt      = min(dx,dy)/1.5;

q   = zeros(N,M);
q(16:32,48:52) = 1;
ax    = ones(N,M)*1;
ay    = zeros(N,M)*0.1;

q(1:48,1:32) = NaN;
ax(1:48,1:32)  = NaN;
ay(1:48,1:32)  = NaN;

qS =  zeros(N,1);
qN =  zeros(N,1);
qW =  zeros(1,M);
qE =  zeros(1,M);


for i=1:50
    
    %...semi-Lagrangian advection
    q    = SemiLagrAdvectMain(ax,ay,q,qS,qN,qW,qE);
    
    
    if (mod(i,1)==0)
        fprintf('time step %i \n',i)
        %...graphics output
        figure(1);
        surf(xx,yy,q,'EdgeColor','None')
        axis image
        axis([0 xLen 0 yLen]);
        drawnow
        %pause(0.1)
    end
end
