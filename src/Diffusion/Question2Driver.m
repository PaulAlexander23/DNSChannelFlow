
%--------------------------------------
% Driver routine for Question two
%--------------------------------------

clear all
close all

global Rd Ad Pd
global xLen yLen
global dx dy
global N M
global D dt

%...diffusion coefficient
D = 1/500;

%...domain size
xLen  = 0.80;
yLen  = 0.64;

%...resolution, fine-mesh spacing
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

%...setup the hierarchy of R,A,P
getRPd(64,32);
getAd(64,32);

%...setup potential field, boundary conditions
q    = zeros(N,M);
q(1:48,1:32)  = NaN;


uS =  zeros(N,1);
uN =  ones(N,1);
uW =  zeros(1,M);
uE =  zeros(1,M);

Qbc      = zeros(N,M);
Qbc(1,33:64) = uW(33:64)*dt*D/dx/dx;
Qbc(49,1:32) = uW(1:32)*dt*D/dx/dx;
Qbc(N,:) = uE*dt*D/dx/dx;
Qbc(1:48,33) = uS(1:48)*dt*D/dy/dy;
Qbc(49:80,1) = uS(49:80)*dt*D/dy/dy;
Qbc(:,M) = uN*dt*D/dy/dy;

for i=1:150
    
    %...diffusion
    q    = DiffusionSolveMain(q,Qbc);
    
    if (mod(i,10)==0)
        fprintf('time step %i \n',i)
        %...graphics output
        figure(1);
        surf(xx,yy,q,'EdgeColor','None')
        colorbar;
        axis image
        axis([0 xLen 0 yLen]);
        title('A plot of the Solution at Time = 1')
        xlabel('x')
        ylabel('y')
        drawnow
    end
end
