A = rand(1,12);
disp(A');
E = cat(2,eye(3), zeros(3,1));
%disp(E);
F = kron(eye(4),E);
%disp(F)
v = A*F;
%disp(v);

B = ones(1,4);
disp(B');
G = zeros(1,4);
G(1,4) = 1;
%disp(G)
H = kron(eye(4),G);
%disp(H)
w = B*H;
%disp(w);

u = v + w;
disp(u')
