function x = CG(AJ,d,x0,tol,scale)
    x = x0;
    AJT = AJ';
    r = d- (AJT*(AJ*x) + scale*x);
    p = r;
    r_square = r'*r;
    while(r_square > tol)
        r_square_before = r_square;
        Hp = AJT*(AJ*p) + scale*p;
        pHp = p'*Hp;
        alpha = r_square_before/pHp;
        x = x + alpha*p;
        r = r - alpha*Hp;
        r_square = r'*r;
        beta = r_square/r_square_before;
        p = r+beta*p;
    end
    
end