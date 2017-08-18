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
  xLen  = 0.8;
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

  

  p0   = zeros(N,M); 
  u    = zeros(N,M); 
  v    = zeros(N,M);
  p0(1:48,1:32) = NaN;
  u(1:48,1:32)  = NaN;
  v(1:48,1:32)  = NaN;
  
  u(30:80,40:50) = 1; 
  
  
  %...setup the hierarchy of R,A,P
  getRPd(64,32);
  getAd(64,32);

  uS =  zeros(N,1);
  uN =  zeros(N,1);
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
  

  for i=1:10

    %...semi-Lagrangian advection
    %uast    = SemiLagrAdvectMain(u,v,u,uS,uN,uW,uE);
    %vast    = SemiLagrAdvectMain(u,v,v,vS,vN,vW,vE);
    
    %...diffusion
    unew    = DiffusionSolveMain(u,Ubc);
    vnew    = DiffusionSolveMain(v,Vbc);
    
    u       = unew;
    v       = vnew;
    
    %...computing divergence
    %[ux,uy] = Diff(unew,N,M,'D',uS,uN,uW,uE);
    %[vx,vy] = Diff(vnew,N,M,'D',vS,vN,vW,vE);
    %Div     = ux + vy; Div = Div(:); Div(1) = 0;

    %...solving for pressure
    %pnew    = PoissonSolve(p0,Div);
    %p0      = pnew;

    %...correcting velocities
    %[px,py] = Diff(pnew,N,M,'N',vS,vN,vW,vE);
    %u       = unew - px(:);
    %v       = vnew - py(:);

    if (mod(i,1)==0)     
      k = 2;
      fprintf('time step %i \n',i)
      %...graphics output
      [ux,uy] = Diff(u,N,M,'D',uS,uN,uW,uE);
      [vx,vy] = Diff(v,N,M,'D',vS,vN,vW,vE);
      vort    = uy-vx;
      vort2D  = reshape(vort,N,M);
      u2D     = reshape(u,N,M);
      v2D     = reshape(v,N,M);
      vel  = sqrt(u2D.*u2D + v2D.*v2D);

      figure(5);
      %quiver(xx(1:k:end,1:k:end),yy(1:k:end,1:k:end),u2D(1:k:end,1:k:end),v2D(1:k:end,1:k:end),12/k)
      contourf(xx,yy,vel,100,'EdgeColor','None')
      %contourf(xx,yy,vort2D,100,'EdgeColor','None')
      colorbar;
      axis image
      axis([0 xLen 0 yLen]);
      drawnow
    end
  end
