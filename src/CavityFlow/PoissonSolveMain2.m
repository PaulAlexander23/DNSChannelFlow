function pnew = PoissonSolveMain2(p,rhs)

    global dx dy
    
    p1 = p(1:64,33:64);
    p2 = p(49:80,1:64);
    error = ones(1,101);
    
    i = 0;
    while (error(i+1)>1e-9 && i < 5)
        B1 = zeros(64,32);
        B1(49:64,1) = -(p2(1:16,32) + p2(1:16,33))/dy/dy;
        B1(64,1:32) = B1(64,1:32) - (p2(16,33:64) + p2(17,33:64))/dx/dx;
        p1 = PoissonSolve(p(1:64,33:64),rhs(1:64,33:64)+B1,1);        
        
        B2 = zeros(32,64);
        B2(1,33:64) = -(p1(48,1:32) + p1(49,1:32))/dx/dx;
        p2 = PoissonSolve(p(49:80,1:64)',(rhs(49:80,1:64)+B2)',2)';
        
        i = i + 1;
        temp = max(max(abs(p1(49:64,:)-p2(1:16,33:64))));
        error(i) = temp;  
    end
    
    disp(i)
    semilogy(error(2:i));
    title('A Plot showing the Decrease in Residual on the Overlap Region');
    xlabel('Iteration');
    ylabel('Maximum Difference on the Overlap Region');
    fprintf('Error: %g \n',error(i))
    fprintf('Rate: %g \n',error(i-1)/error(i))
        

    pnew = zeros(80,64);
    pnew(1:64,33:64) = p1;
    pnew(49:80,1:64) = p2;
    pnew(1:48,1:32) = NaN;
    
end

