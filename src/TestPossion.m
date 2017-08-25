
%--------------------------------------
% Driver routine for Question three
%--------------------------------------

clear all
close all

global Rp Ap Pp
global xLen yLen
global dx dy
global N M
global dt

%...domain size
xLen  = 2;
yLen  = 1;

%...resolution, fine-mesh spacing
N     = 64;
M     = 32;

x       = linspace(0,xLen,N+1);
y       = linspace(0,yLen,M+1);
xc      = (x(1:end-1)+x(2:end))/2;
yc      = (y(1:end-1)+y(2:end))/2;
[yy,xx] = meshgrid(yc,xc);
dx      = xLen/N;
dy      = yLen/M;
dt      = min(dx,dy)/1.5;

p0   = -1.4*ones(N,M);

%...setup the hierarchy of R,A,P
getRPp(64,32);
getAp(64,32);

f = @(x,y) 1 + 0.*x;
Div = f(xx,yy);


p1 = PoissonSolve(p0,Div);

%...graphics output
figure;
surf(xx,yy,p1,'EdgeColor','None')
colorbar;
axis image
axis([0 xLen 0 yLen]);
title('The Numerical Solution')
xlabel('x')
ylabel('y')

