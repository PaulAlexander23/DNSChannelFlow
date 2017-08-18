function [qx,qy] = Diff( q,N,M,qS,qN,qW,qE,typ)

global dx dy

qG = zeros(N+2,M+2);
qG(2:N+1,2:M+1)   =  q;

if (nargin == 8)
    qG(1,34:65)    = q(1,33:64) + 2*typ(1,33:64).*(qW(33:64)-q(1,33:64));
    qG(49,2:33)    = q(49,1:32) + 2*typ(1,1:32).*(qW(1:32)-q(49,1:32));
    qG(82,2:65)    = q(80,:)    + 2*typ(80,:).*(qE-q(80,:));
    qG(2:49,33)    = q(1:48,33) + 2*typ(1:48,33).*(qS(1:48)-q(1:48,33));
    qG(50:81,1)    = q(49:80,1) + 2*typ(49:80,1).*(qS(49:80)-q(49:80,1));
    qG(2:81,66)    = q(:,64)    + 2*typ(:,64).*(qN-q(:,64));
else
    qG(1,34:65)   = 2*qW(33:64)-q(1,33:64);
    qG(49,2:33)   = 2*qW(1:32)-q(49,1:32);
    qG(82,2:65)   = 2*qE-q(80,:);
    qG(2:49,33)   = 2*qS(1:48)-q(1:48,33);
    qG(50:81,1)   = 2*qS(49:80)-q(49:80,1);
    qG(2:81,66)   = 2*qN-q(:,64);
end

qx = zeros(N,M);
qy = zeros(N,M);
qx(1:80,33:64) = (qG(3:82,34:65)-qG(1:80,34:65))/(2*dx);
qx(49:80,1:32) = (qG(51:82,2:33)-qG(49:80,2:33))/(2*dx);
qy(1:48,33:64) = (qG(2:49,32+3:M+2)-qG(2:49,32+1:M))/(2*dy);
qy(49:80,1:64) = (qG(50:81,3:M+2)-qG(50:81,1:M))/(2*dy);
%qx = (qG(3:N+2,2:M+1)-qG(1:N,2:M+1))/(2*dx);
%qy = (qG(2:N+1,3:M+2)-qG(2:N+1,1:M))/(2*dy);
qx(1:48,1:32)  = NaN;
qy(1:48,1:32)  = NaN;


end

