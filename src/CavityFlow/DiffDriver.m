N=80;
M=64;
u    = zeros(N,M);
%u(1:48,33:64)  = 1;
u = -(xx.^2 + yy.^2);
u(1:48,1:32)  = NaN;

uS =  zeros(N,1);
uN =  zeros(N,1);
uE =  ones(1,M);
uW =  zeros(1,M);

typ = zeros(80,64);

[ux,uy] = Diff(u,N,M,uS,uN,uW,uE,typ);

figure;
surf(xx,yy,u,'EdgeColor','None');
%axis image
axis([0 xLen 0 yLen]);
figure;
surf(xx,yy,ux,'EdgeColor','None');
%axis image
axis([0 xLen 0 yLen]);
figure;
surf(xx,yy,uy,'EdgeColor','None');
%axis image
axis([0 xLen 0 yLen]);