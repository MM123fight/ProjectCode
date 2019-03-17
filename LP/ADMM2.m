%read data
root = "./data";
experiment_names = {"Mat0|1000_10000_0.1_0.5"};
experiment = root + "/" + experiment_names{1};
[nb,nf,mi,me,m,n,Ae,Ai,be,bi,c] = dataRead(experiment);

%parameter
rho = 1;
tol = 1e-15;
sub_tol = 1e-3;
grad_tol = 2*sub_tol*rho;
mu = 0.49;
delta = 0.9;

%ADMM1
m_bar = 2*me +mi;
A = [Ai;Ae;-Ae];
b = [bi;be;-be];
AT = A';
H = A*AT + eye(m_bar);
y = zeros(n,1);
x = y;
z = zeros(m_bar+n,1);
condition = 0.1;
tic
while( condition > tol)
    v = x + (c+z(m_bar+1:m_bar+n))/rho - y;
    w = A*x - b + z(1:m_bar)/rho;
    J = find(w>0);
    wJ = w(J);
    AJ = A(J,:);
    AJT = AJ';
    grad = rho*(AJT*wJ +v);
    grad_cond = grad'*grad;
    %SSN however, not always use semismooth
    while(grad_cond > grad_tol)
        h = AJ*v - w(J);
        dir = AJT*(H(J,J)\h) - v;
        Adir = A*dir;
        tmp = -0.5*rho*(wJ'*wJ);
        tmp1 = dir'*(c + rho*v- mu*grad);
        tmp2 = 0.5*rho*(dir'*dir);
        eta_before = 0;
        %line search
        for i = 0:20
            eta = power(delta,i);
            w = w + (eta-eta_before)*Adir;
            J = find(w > 0);
            wJ = w(J);
            line_search = 0.5*rho*(wJ'*wJ) + tmp+ eta*tmp1 + power(eta,2)*tmp2;
            if(line_search <= 0)
                break;
            end
            eta_before = eta;
        end
        v = v + eta*dir;
        x = x + eta*dir;
        AJ = A(J,:);
        AJT = AJ';
        grad = rho*(AJT*wJ +v);
        grad_cond = grad'*grad;
    end
    y = x + z(m_bar+1:m_bar+n)/rho;
    y(1:nb) = max(y(1:nb),0);
    d = A*x - b;
    xminusy = x - y;
    z(1:m_bar) = max(z(1:m_bar)+rho*d,0);
    z(1+m_bar:m_bar+n) = z(1+m_bar:m_bar+n) + rho*xminusy;
    s = max(d,0);
    condition = max(max(s.*s),max(xminusy.*xminusy)); 
end
toc