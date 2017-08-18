%--------------------------------------
% Driver routine for cavity flow
%--------------------------------------

  clear all
  close all

  global xLen
  global dx
  global N
  global dt Re

  %...Reynolds number
  Re = 500;

  %...domain size
  xLen  = 2*pi;
  
  N = 150;
  
  x       = linspace(0,xLen,N+1)';
  xc      = (x(1:end-1)+x(2:end))/2;
  dx      = xLen/N;
  dt      = 0.01;

  u    = 0.1*ones(N,1); 
  h    = cos(xc)*0.1+0.5;

  %u = 0.05*ones(N,1); 
  %h = 0.5*ones(N,1); 
  %h(1)=0.6;
  
  f_h = @(u,h,du,dh) u.*dh+h.*du;
  f_u = @(u,h,du,dh) (1-3*h)./(1-h).*u.*du + ((1-h)-1./(1-h.^2).*u.^2).*dh;
  

   figure(1);
      plot(xc,[h,u]);
      axis([0 xLen 0 1]);
      drawnow
      pause(1);
      
  for i=1:6000

    %...semi-Lagrangian advection
    u    = Advect(u,u,h,f_u);
    h    = Advect(h,u,h,f_h);
    

    if (mod(i,10)==0)
      fprintf('time step %i \n',i)
      %...graphics output
      figure(1);
      plot(xc,[h,u,u.^2./(1-h).^2/2,0.5*ones(N,1)]);
      axis([0 xLen 0 1]);
      drawnow
      pause(0.01);
    end
  end
