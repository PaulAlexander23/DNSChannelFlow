%--------------------------------------
% Driver routine for cavity flow
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
  xLen  = 80;
  yLen  = 32;

  %...resolution, fine-mesh spacing
  %   and number of cycles
  p     = 6;
  q     = 6;
  ncyc  = 10;

  N     = 80;
  M     = 32;

  x       = linspace(0,xLen,N+1);
  y       = linspace(0,yLen,M+1);
  xc      = (x(1:end-1)+x(2:end))/2;
  yc      = (y(1:end-1)+y(2:end))/2;
  [yy,xx] = meshgrid(yc,xc);
  dx      = xLen/N;
  dy      = yLen/M;
  dt      = min(dx,dy)/1.5;

  %...setup the hierarchy of R,A,P
%   getRPd(N,M);
%   getRPp(N,M);
%   getAd(N,M);
%   getAp(N,M);

  p0   = zeros(N,M); 
  u    = ones(N,M)*0.1; 
  v    = zeros(N,M);
  
  u(10:20,10:20) = 0.1;
  v(10:20,10:20) = 0;
  
  p0 = p0(:);
  u  = u(:);
  v  = v(:);

  uS =  ones(N,1)*0.1;
  uN =  ones(N,1)*0.1;
  uW =  ones(1,M)*0.1;
  uE =  ones(1,M)*0.1;

  vS = zeros(N,1);
  vN = zeros(N,1);
  vW = zeros(1,M);
  vE = zeros(1,M);

  Ubc      = zeros(N,M);
  Ubc(1,:) = uW*dt/Re/dx/dx;
  Ubc(N,:) = uE*dt/Re/dx/dx;
  Ubc(:,1) = uS*dt/Re/dy/dy;
  Ubc(:,M) = uN*dt/Re/dy/dy;
  Ubc      = Ubc(:);

  Vbc      = zeros(N,M);
  Vbc(1,:) = vW*dt/Re/dx/dx;
  Vbc(N,:) = vE*dt/Re/dx/dx;
  Vbc(:,1) = vS*dt/Re/dy/dy;
  Vbc(:,M) = vN*dt/Re/dy/dy;
  Vbc      = Vbc(:);

  for i=1:1500

    %...semi-Lagrangian advection
    uast    = SemiLagrAdvect(u,v,u,uS,uN,uW,uE);
    vast    = SemiLagrAdvect(u,v,v,vS,vN,vW,vE);
    
    u = uast;
    v = vast;
    

    if (mod(i,10)==0)
      k = 2;
      fprintf('time step %i \n',i)
      %...graphics output
      u2D     = reshape(u,N,M);
      v2D     = reshape(v,N,M);

      figure(1);
      quiver(xx(1:k:end,1:k:end),yy(1:k:end,1:k:end),u2D(1:k:end,1:k:end),v2D(1:k:end,1:k:end),12/k)
      axis image
      axis([0 xLen 0 yLen]);
      drawnow
    end
  end
