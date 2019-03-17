function x = CG_inv(size,AJ,d,tol,scale)
    inv_tol = scale^2*tol;
    z = zeros(size,1);
    AJd = AJ*d;
    AJT = AJ';
    r = AJd- (AJ*(AJT*z) + scale*z);
    p = r;
    r_square = r'*r;
    AJTp = AJT*p;
    res = AJTp;
    res_square = res'*res;
    while(res_square > inv_tol)
        r_square_before = r_square;
        Hp = AJ*AJTp + scale*p;
        pHp = p'*Hp;
        alpha = r_square_before/pHp;
        z = z + alpha*p;
        r = r - alpha*Hp;
        r_square = r'*r;
        beta = r_square/r_square_before;
        res = AJT*r;
        p = r+beta*p;
        AJTp = res+beta*AJTp;
        res_square = res'*res;
    end
    x = scale*(d-AJT*z);
    
end
