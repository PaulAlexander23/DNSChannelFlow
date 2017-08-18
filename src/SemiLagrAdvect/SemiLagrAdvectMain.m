function qnew = SemiLagrAdvectMain(u,v,q,qS,qN,qW,qE)

    q1 = SemiLagrAdvectModified(u(1:48,33:64),v(1:48,33:64),q(1:48,33:64),qS(1:48),qN(1:48),qW(33:64),(q(48,33:64)+q(49,33:64))*0.5,48,32);
    q2 = SemiLagrAdvectModified(u(49:80,1:64),v(49:80,1:64),q(49:80,1:64),qS(49:80),qN(49:80),[qW(1:32),(q(48,33:64)+q(49,33:64))*0.5],qE,32,64);
    
    qnew = zeros(80,64);
    qnew(1:48,33:64) = q1;
    qnew(49:80,1:64) = q2;
    qnew(1:48,1:32) = NaN;
    
end
