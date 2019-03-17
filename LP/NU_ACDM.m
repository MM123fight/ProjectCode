function x = NU_ACDM(n,A,d,scale,x0,tol,cond_type)
    %(AT*A+scale*I)x = d
    %parameter
    beta = 0;
    sigma_beta = scale;
    alpha = (1-beta)/2;
    L = sum(A.*A,1)' + scale*ones(n,1);
    Lalpha = power(L,alpha);
    Salpha = sum(Lalpha);
    p = Lalpha/Salpha;
    tau = 2/(1+sqrt(4*Salpha^2/sigma_beta +1));
    eta = 1/(tau*Salpha^2);
    eta_sigma = eta*sigma_beta;
    
    %variables
    y = x0;
    z = x0;
    x = x0;
    w = A*x;
    G = A'*w + x - d;
    G_square = G.*G;
    if(cond_type == "grad_cord")
        cond = max(G_square,tol);
    else 
        cond = G_square;
    end
    sum_cond = sum(cond);
    while( sum_cond > n*tol)
       %sample i with prob ~ (sqrt(Li))
       i = randsrc(1,1,[1:n; p']);

       grad = A(:,i)'*w +scale*x(i)- d(i);
       y(i) = x(i) - 1/L(i)*grad;
       z(i) = 1/(1+eta_sigma)*(z(i)+eta_sigma*x(i)-eta/p(i)*grad);
       tmp = tau*z(i) + (1-tau) *y(i);
       delta_x = tmp - x(i);
       x(i) = tmp;
       w = w + delta_x*A(:,i);
  
       cond_tmp = grad^2;
       if(cond_type == "grad_cord")
           cond_tmp = max(grad^2,tol);
       end
       sum_cond = sum_cond + cond_tmp - cond(i);
       cond(i) = cond_tmp;
    end
end