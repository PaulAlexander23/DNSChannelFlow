%------------------------------------------------
% Test routine for advection across the boundary
%------------------------------------------------

  clear all
  close all

  addpath ~/Desktop/FOR_APS/NumPDE_Codes/colormaps/

  global Rd Ad Pd
  global Rp Ap Pp
  global xLen yLen
  global dx dy
  global N M
  global dt Re

  %...Reynolds number
  Re = 500;

  %...domain size
  xLen  = 2;
  yLen  = 1;

  %...resolution, fine-mesh spacing
  %   and number of cycles
  p     = 7;
  q     = 6;
  ncyc  = 10;

  N     = 2^p;
  M     = 2^q;

  x       = linspace(0,xLen,N+1);
  y       = linspace(0,yLen,M+1);
  xc      = (x(1:end-1)+x(2:end))/2;
  yc      = (y(1:end-1)+y(2:end))/2;
  [yy,xx] = meshgrid(yc,xc);
  dx      = xLen/N;
  dy      = yLen/M;
  dt      = min(dx,dy)/1.5;

  %...setup the hierarchy of R,A,P
  getRPd(N,M);
  getRPp(N,M);
  getAd(N,M);
  getAp(N,M);

  p0   = zeros(N,M); p0 = p0(:);
  u    = 0.2*ones(N,M); u(1:N*0.25,M/4:3*M*0.25) = 1; u = u(:);
  v    = zeros(N,M); v = v(:);

  uS =  zeros(N,1);
  uN =  zeros(N,1);
  uW =  zeros(1,M);
  uE =  zeros(1,M);

  vS = zeros(N,1);
  vN = zeros(N,1);
  vW = zeros(1,M);
  vE = zeros(1,M);

  Ubc      = zeros(N,M);
  Ubc(1,:) = uW*dt/2/Re/dx/dx;
  Ubc(N,:) = uE*dt/2/Re/dx/dx;
  Ubc(:,1) = uS*dt/2/Re/dy/dy;
  Ubc(:,M) = uN*dt/2/Re/dy/dy;
  Ubc = Ubc(:);

  Vbc      = zeros(N,M);
  Vbc(1,:) = vW*dt/2/Re/dx/dx;
  Vbc(N,:) = vE*dt/2/Re/dx/dx;
  Vbc(:,1) = vS*dt/2/Re/dy/dy;
  Vbc(:,M) = vN*dt/2/Re/dy/dy;
  Vbc = Vbc(:);

  for i=1:1000

    %...semi-Lagrangian advection
    u    = SemiLagrAdvect(u,v,u,uS,uN,uW,uE);
    v    = SemiLagrAdvect(u,v,v,vS,vN,vW,vE);

    if (mod(i,20)==0)
      k = 8;
      fprintf('time step %i \n',i)
      %...graphics output
      [ux,uy] = Diff(u,N,M,'D',uS,uN,uW,uE);
      [vx,vy] = Diff(v,N,M,'D',vS,vN,vW,vE);
      vort    = uy-vx;
      vort2D  = reshape(vort,N,M);
      u2D  = reshape(u,N,M);
      v2D  = reshape(v,N,M);
      vel  = sqrt(u2D.*u2D + v2D.*v2D);

      figure(1);
%      quiver(xx(1:k:end,1:k:end),yy(1:k:end,1:k:end),u2D(1:k:end,1:k:end),v2D(1:k:end,1:k:end),15/k)
      contourf(xx,yy,vel,100,'EdgeColor','None')
      axis image
      axis([0 xLen 0 yLen]);
      colorbar
      drawnow
    end
  end
