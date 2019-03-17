%read data
root = "./data";
experiment_names = {"Mat0|100_1000_0.1_0.5"};
experiment = root + "/" + experiment_names{1};
[nb,nf,mi,me,m,n,Ae,Ai,be,bi,c] = dataRead(experiment);



%parameter for outer iterations
rho = 1;
tol = 1e-3;
tau = 0.5*(1 + sqrt(5));
sigma = tau*rho;
%parameter for inner iterations
sub_tol = 1e-3;
ACDMe_cond = "grad_cord";
ACDMe_tol = 2*sub_tol;
SSN_tol = 2*sub_tol;
mu = 0.49;
delta = 0.9;


%primal variables
x = zeros(n,1);
y = zeros(n,1);
z = zeros(n,1);
%dual variables
lambda_I = zeros(mi,1);
lambda_E = zeros(me,1);
lambda_x = zeros(n,1);
lambda_z = zeros(n,1);

exact_e = check_me(me);

AeT = Ae';
if(exact_e == true)
    Inv =  inv(Ae*AeT+ eye(me));
end
%ADMM
condition = 1;
tic;
while(condition > tol)
    %step1:update x and z
    if(exact_e == true)
        %update x: exactly
        we = lambda_E/rho - be;
        ve = (lambda_x+c)/rho - y;
        x = -ve +AeT*(Inv*(Ae*ve-we));
        fprintf("check1");
    else
        %inexactly -> ACDM ;
        h = -(AeT*we+ve);
        x = NU_ACDM(n,Ae,h,1,x,ACDMe_tol,ACDMe_cond);
        fprintf("check2");
    end
    
    %update z: inexactly -> SSN(Linear System:ACDM?conjugate gradient?)
    wi = Ai*z+lambda_I/rho-bi;
    vi = z+(lambda_z +c)/rho-y;
    fprintf("check3");
    z = SSN(mi,n,Ai,z,wi,vi,SSN_tol,c,rho);
    fprintf("check4");
    
    %step2:update y
    y = 0.5*(x+z+(lambda_x + lambda_z)/rho);
    y(1:nb) = max(y(1:nb),0);
    
    %step3:update lambda
    de = Ae*x - be;
    di = Ai*z - bi;
    xminusy = x-y;
    zminusy = z-y;
    lambda_E = lambda_E + sigma*de;
    lambda_I = (1-tau)*lambda_I + tau*max( rho* di + lambda_I,0);
    lambda_x = lambda_x + sigma*xminusy;
    lambda_z = lambda_z + sigma*zminusy;
     
    %check the stopping condition
    gap = c'*x + be'*lambda_E + bi'*lambda_I;
    tmp = max(di,0); 
    condition = max(max(de.*de),max(tmp.*tmp));
    condition = max(condition,max(xminusy.*xminusy));
    condition = max(condition,max(zminusy.*zminusy));
    condition
end
toc

