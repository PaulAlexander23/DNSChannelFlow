function u = MGVp(ilev,u,rhs,i)
%======================================
% multigrid V-cycle (recursive)
%======================================

  global Rp Ap Pp
  global ilevmin

  mu = 5;

  %...lowest level
  if (ilev==ilevmin+1)
    u  = GSp(ilev,u,rhs,mu,i);
  else
    %...presmoothing
    u  = GSp(ilev,u,rhs,mu,i);

    %...residual and restriction
    r  = rhs - Ap{ilev,i}*u;
    rr = Rp{ilev,i}*r;

    %...call MGV
    ee = zeros(size(rr));
    ee = MGVp(ilev+1,ee,rr,i);

    %...prolong and add
    u  = u + Pp{ilev,i}*ee;

    %...postsmoothing
    u  = GSp(ilev,u,rhs,mu,i);
  end
