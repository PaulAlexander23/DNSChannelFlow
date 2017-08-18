function qnew = DiffusionSolveMain(q,Qbc)

global dx dy
global dt D

q1 = q(1:64,33:64);
q2 = q(49:80,1:64);
error = 1;
rate = 2;
i = 0;

while error>1e-6 && rate>1.1 && i < 10
    
    B1 = zeros(64,32);
    B1(49:64,1) = (q2(1:16,32) + q2(1:16,33))*dt*D/dy/dy/2;
    B1(64,1:32) = B1(64,1:32) + (q2(16,33:64) + q2(17,33:64))*dt*D/dx/dx/2;
    q1 = DiffusionSolveModified(q(1:64,33:64),Qbc(1:64,33:64)+B1);
    
    B2 = zeros(32,64);
    B2(1,33:64) = (q1(48,1:32) + q1(49,1:32))*dt*D/dx/dx/2;
    q2 = DiffusionSolveModified(q(49:80,1:64)',(Qbc(49:80,1:64)+B2)')';
    
    temp = max(abs(q1(49:64,1)-q2(1:16,33)));
    temp = max(temp,max(abs(q1(64,:)-q2(16,33:64))));
    temp = max(temp,max(abs(q1(49,:)-q2(1,33:64))));
    rate = error/temp;
    error = temp;
    
    i = i + 1;
end

qnew = zeros(80,64);
qnew(1:64,33:64) = q1;
qnew(49:80,1:64) = q2;
qnew(1:48,1:32) = NaN;

end

