function qnew = DiffusionSolveMain2(q,Qbc)

    global dx dy
    global dt D

    q1 = q(1:64,33:64);
    q2 = q(49:80,1:64);
    error = ones(1,11);
    rate = 2*ones(1,11);
    i = 2;
    while error(:,i-1)>1e-10 && max(rate(:,i-1))>1.1 && i < 11

        B1 = zeros(64,32);
        B1(49:64,1) = (q2(1:16,32) + q2(1:16,33))*dt*D/dy/dy/2;
        B1(64,1:32) = B1(64,1:32) + (q2(16,33:64) + q2(17,33:64))*dt*D/dx/dx/2;
        q1 = DiffusionSolveModified(q(1:64,33:64),Qbc(1:64,33:64)+B1);
        
        B2 = zeros(32,64);
        B2(1,33:64) = (q1(48,1:32) + q1(49,1:32))*dt*D/dx/dx/2;
        q2 = DiffusionSolveModified(q(49:80,1:64)',(Qbc(49:80,1:64)+B2)')';

        
        temp = max(max(abs(q1(49:64,:)-q2(1:16,33:64))));
        rate(i) = error(i-1)/temp;
        error(i) = temp;        
        
        fprintf(1,'norm1: %g ',error(:,i));
        fprintf(1,'rate1: %g \n',rate(:,i));
        i = i + 1;
    end
    
    figure(2);
    semilogy(error(2:i-1));
    title('A Plot Demonstrating the Decrease of the Residual on the Interface')
    xlabel('Iteration')
    ylabel('Maximum Difference on the Overlap Region')
    
%     figure(3);
%     plot(max(rate(:,3:i)));
    
    disp(i)
    qnew = zeros(80,64);
    qnew(1:64,33:64) = q1;
    qnew(49:80,1:64) = q2;
    qnew(1:48,1:32) = NaN;
    
end

