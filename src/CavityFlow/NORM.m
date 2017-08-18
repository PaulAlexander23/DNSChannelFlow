function rn = NORM(p0,rhs,i)
    global Ap


    p0 = p0(:);
    rhs = rhs(:);

    rn  = norm(rhs - Ap{1,i}*p0);

end

