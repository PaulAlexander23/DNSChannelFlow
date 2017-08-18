function qnew = SemiLagrAdvectModified(u,v,q,qS,qn,qW,qE,n,m)

global dx dy
global dt Re

%...embedding
qq              = zeros(n+2,m+2);
qq(2:n+1,2:m+1) = q;

%...set the ghost values (four edges)
qq(1,2:m+1)   = 2*qW-qq(2,2:m+1);
qq(n+2,2:m+1) = 2*qE-qq(n+1,2:m+1);
qq(2:n+1,1)   = 2*qS-qq(2:n+1,2);
qq(2:n+1,m+2) = 2*qn-qq(2:n+1,m+1);

%...set the ghost values (four corners)
qq(1,1)       = -qq(2,2);
qq(n+2,1)     = -qq(n+1,2);
qq(n+2,m+2)   = -qq(n+1,m+1);
qq(1,m+2)     = -qq(2,m+1);

q1   = qq(2:n+1,2:m+1);

q2p  = qq(3:n+2,2:m+1);
q2m  = qq(1:n,2:m+1);

q3p  = qq(2:n+1,3:m+2);
q3m  = qq(2:n+1,1:m);

q4pp = qq(3:n+2,3:m+2);
q4mm = qq(1:n,1:m);
q4pm = qq(3:n+2,1:m);
q4mp = qq(1:n,3:m+2);

xi  = -u*dt/dx;
eta = -v*dt/dy;

Q2  = q2p.*(xi>0) + q2m.*(xi<0);
Q3  = q3p.*(eta>0) + q3m.*(eta<0);
Q4  = q4pp.*((xi>0) & (eta>0)) + q4mm.*((xi<0) & (eta<0)) + ...
    q4pm.*((xi>0) & (eta<0)) + q4mp.*((xi<0) & (eta>0));

qnew = (1-abs(xi)).*(1-abs(eta)).*q1 + ...
    abs(xi).*(1-abs(eta)).*Q2 + ...
    abs(eta).*(1-abs(xi)).*Q3 + ...
    abs(xi).*abs(eta).*Q4;


