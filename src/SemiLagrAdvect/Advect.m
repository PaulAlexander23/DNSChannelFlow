function qnew = Advect( q,u,h,f)


global N
global dx
global dt

%...embedding
qq          = zeros(N+2,1);
qq(2:N+1)   = q;
uu          = zeros(N+2,1);
uu(2:N+1)   = u;
hh          = zeros(N+2,1);
hh(2:N+1)   = h;

%...set the ghost values (four corners)
uu(1)   = uu(N+1);
uu(N+2) = uu(2);
hh(1)   = hh(N+1);
hh(N+2) = hh(2);

du = (uu(3:N+2)-uu(1:N))/(2*dx);
dh = (hh(3:N+2)-hh(1:N))/(2*dx);

qnew = q + dt*f(u,h,du,dh);

end

