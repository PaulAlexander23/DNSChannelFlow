%--------------------------------------
% Driver routine for Question four
%--------------------------------------

%clear all
close all

global Rd Ad Pd
global Rp Ap Pp
global xLen yLen
global dx dy
global N M
global dt Re

%...Reynolds number
Re = 2000;

%...domain size
xLen  = 0.80;
yLen  = 0.64;

%...resolution, fine-mesh spacing
%   and number of cycles
ncyc  = 10;
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
getRPp(64,32);
getAd(64,32);
getAp(64,32);

p0   = zeros(N,M); p0(1:48,1:32) = NaN;
u    = zeros(N,M); u(1:48,1:32)  = NaN;
v    = zeros(N,M); v(1:48,1:32)  = NaN;

uS =  zeros(N,1);
uN =  ones(N,1);
uW =  zeros(1,M);
uE =  zeros(1,M);

vS = zeros(N,1);
vN = zeros(N,1);
vW = zeros(1,M);
vE = zeros(1,M);

Ubc      = zeros(N,M);
Ubc(1,33:64) = uW(33:64)*dt/Re/dx/dx;
Ubc(49,1:32) = uW(1:32)*dt/Re/dx/dx;
Ubc(N,:) = uE*dt/Re/dx/dx;
Ubc(1:48,33) = uS(1:48)*dt/Re/dy/dy;
Ubc(49:80,1) = uS(49:80)*dt/Re/dy/dy;
Ubc(:,M) = uN*dt/Re/dy/dy;

Vbc      = zeros(N,M);
Vbc(1,33:64) = vW(33:64)*dt/Re/dx/dx;
Vbc(49,1:32) = vW(1:32)*dt/Re/dx/dx;
Vbc(N,:) = vE*dt/Re/dx/dx;
Vbc(1:48,33) = vS(1:48)*dt/Re/dy/dy;
Vbc(49:80,1) = vS(49:80)*dt/Re/dy/dy;
Vbc(:,M) = vN*dt/Re/dy/dy;


% figure(1);
% quiver(xx(1:2:end,1:2:end),yy(1:2:end,1:2:end),u(1:2:end,1:2:end)0,...
%  v(1:2:end,1:2:end),6)
% axis image
% axis([0 xLen 0 yLen]);
% drawnow

error = ones(1,5001);
rate = 2*ones(1,5000);

for i=1:5000
    %clc;
    %fprintf('time step %i \n',i)
    
    %...semi-Lagrangian advection
    uast    = SemiLagrAdvectMain(u,v,u,uS,uN,uW,uE);
    vast    = SemiLagrAdvectMain(u,v,v,vS,vN,vW,vE);
    
    %...diffusion
    unew    = DiffusionSolveMain(uast,Ubc);
    vnew    = DiffusionSolveMain(vast,Vbc);
    
    %...computing divergence
    %     [ux,uy] = Diff(unew,N,M,uS,uN,uW,uE);
    %     [vx,vy] = Diff(vnew,N,M,vS,vN,vW,vE);
    [ux,uy] = Diff2(unew,N,M,'D',uS,uN,uW,uE);
    [vx,vy] = Diff2(vnew,N,M,'D',vS,vN,vW,vE);
    Div     = ux + vy; Div(1,33) = 0;
    
    %...solving for pressure
    pnew    = PoissonSolveMain(p0,Div);
    p0      = pnew;
    
    %...correcting velocities
    %     [px,py] = Diff(pnew,N,M,vS,vN,vW,vE,typ);
    [px,py] = Diff2(pnew,N,M,'N',vS,vN,vW,vE);
    norm1 = max(max(abs(unew - px - u)));
    norm2 = max(max(abs(vnew - py - v)));
    u       = unew - px;
    v       = vnew - py;
    
    error(i+1) = max(norm1,norm2);
    rate(i) = error(i)/error(i+1);
    
    %fprintf(1,'MaxDiffU: %g MaxDiffV: %g \n',norm1,norm2);
    
    if (mod(i,100)==0)
        k = 2;
        fprintf('time step %i \n',i)
        fprintf(1,'MaxDiffU: %g MaxDiffV: %g \n',norm1,norm2);
        %...graphics output
        figure(1);
        [ux,uy] = Diff2(u,N,M,'D',uS,uN,uW,uE);
        [vx,vy] = Diff2(v,N,M,'D',vS,vN,vW,vE);
        vort    = uy-vx;
        surf(xx,yy,vort,'EdgeColor','None');
        axis image
        axis([0 xLen 0 yLen]);
        colorbar;
        drawnow
        
        figure(2);
        %surf(xx,yy,p0,'EdgeColor','None');
        quiver(xx(1:k:end,1:k:end),yy(1:k:end,1:k:end),u(1:k:end,1:k:end)...
            ,v(1:k:end,1:k:end),12/k)
        axis image
        axis([0 xLen 0 yLen]);
        drawnow
    end
    
end

figure(3);
semilogy(error);
xlabel('Time step')
ylabel('Max Change in Solution')
title('A plot showing the Evolution to the Steady State Solution')

