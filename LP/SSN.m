%function x = SSN(m,n,A,x0,w,v,tol,c,subtype)
    %parameter;
root = "./data";
experiment_names = {"test3","test2","Mat0|10_100_0.1_0.5"};
experiment = root + "/" + experiment_names{1};
[nb,nf,mi,me,m,n,Ae,Ai,be,bi,c] = dataRead(experiment);
mi = mi+nb;
Ai = [Ai;-eye(nb) zeros(nb,n-nb)];
bi = [bi;zeros(nb,1)];
nb = 0; 
m = mi;
A = Ai;
b = bi;

rho = 1;
scale_w = rho;
scale_v = 1/rho;
scale = scale_v/scale_w;
SSN_tau = 0.1;
pow_idx = SSN_tau +1;
sigma = 0.5;
mu = 0.49;
delta = 0.5;
grad_tol = 1e-3;
subtype = "exact";

x0 = zeros(n,1);
x = x0;
%lambda_I = zeros(mi,1);
%lambda_E = zeros(me,1);
lambda = zeros(m,1);
dir = zeros(n,1);
tic 
for j = 1:24
    w = A*x -b + lambda/rho;
    v = zeros(n,1);
    J = find(w>0);
    wJ = w(J);
    AJ = A(J,:);
    AJT = AJ';
    grad =scale_w * AJT*wJ +scale_v*v + c;
    grad_cond = grad'*grad;
    while(grad_cond > grad_tol)
        sizeJ = size(J);
        d = -grad/scale_w;
        if(sizeJ == 0)
            dir = d/scale;
        elseif(subtype == "exact")
            dir = (d - AJT*((AJ*AJT+scale*eye(sizeJ))\(AJ*d)))/scale;
        elseif(subtype == "CG")
            CG_tol = min(sigma,sum(power(d,pow_idx)));
            dir = CG(AJ,d,dir,CG_tol,scale);
        end
        Adir = A*dir;
        tmp = -0.5*scale_w*(wJ'*wJ);
        tmp1 = dir'*(c + scale_v*v+ mu*d);
        tmp2 = 0.5*scale_v*(dir'*dir);
        eta_before = 0;
        %line search
        for i = 0:20
            eta = power(delta,i);
            w = w + (eta-eta_before)*Adir;
            J = find(w > 0);
            wJ = w(J);
            line_search = 0.5*scale_w*(wJ'*wJ) + tmp+ eta*tmp1 + power(eta,2)*tmp2;
            if(line_search <= 0)
                break;
            end
            eta_before = eta;
        end
        v = v + eta*dir;
        x = x + eta*dir;
        AJ = A(J,:);
        AJT = AJ';
        grad = scale_w * AJT*wJ +scale_v*v + c;
        grad_cond = grad'*grad;
        
    end
    lambda = zeros(m,1);
    lambda(J) = rho*wJ;
end
toc
%end