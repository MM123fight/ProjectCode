%read data
root = "./data";
experiment_names = {"test3","Mat10|100_1000_0.1_0.5","Mat0|1000_10000_0.1_0.5"};
experiment = root + "/" + experiment_names{1};
[nb,nf,mi,me,m,n,Ae,Ai,be,bi,c] = dataRead(experiment);

%parameter for outer iterations
rho = 1;
tol = 1e-10;
%tau = 0.5*(1 + sqrt(5));
tau = 1;
sigma = tau*rho;
%subtype = "NU_ACDM";
subtype = "exact";
ACDM_tol = tol;
cond_type = "grad_cord";


%primal variables
x = zeros(n,1);
y = zeros(n,1);
z = zeros(n,1);
Xi = zeros(mi,1);
%dual variables
lambda_I = zeros(mi,1);
lambda_E = zeros(me,1);
lambda_x = zeros(n,1);
lambda_z = zeros(n,1);

time = 0;
tic;
AeT = Ae';
AiT = Ai';
if(subtype == "exact")
    t1 = clock;
    Inve =  inv(Ae*AeT+ eye(me));
    Invi =  inv(Ai*AiT+ eye(mi));
    t2 = clock;
    time = etime(t2,t1);
end
%ADMM
condition = 1;
scale = 1;
iter = 1;
while(condition > tol)
    t1 = clock;
    we = lambda_E/rho - be;
    ve = (lambda_x+c)/rho - y; 
    wi = lambda_I/rho -bi + Xi;
    vi = (lambda_z+c)/rho - y;
    if(subtype == "exact")
        x = -ve +AeT*(Inve*(Ae*ve-we));
        z =  -vi +AiT*(Invi*(Ai*vi-wi));
    elseif(subtype == "NU_ACDM")
        ge = -(AeT*we+ve);
        x = NU_ACDM(n,Ae,ge,scale,zeros(n,1),ACDM_tol,cond_type);
    end
  
    %step2:update y and Xi
    y = 0.5*(x+z+(lambda_x + lambda_z)/rho);
    y(1:nb) = max(y(1:nb),0);
    de = Ae*x - be;
    di = Ai*z - bi;
    tmp_vec = di + lambda_I/rho;
    Xi = -min(tmp_vec,0);
    
    %step3:update lambda
    xminusy = x-y;
    zminusy = z-y;
    lambda_E = lambda_E + sigma*de;
    tmp = di+Xi;
    lambda_I = lambda_I + sigma*tmp;
    lambda_x = lambda_x + sigma*xminusy;
    lambda_z = lambda_z + sigma*zminusy;
    t2 = clock;
    time = time + etime(t2,t1);
     
    %check the stopping condition
    gap =0.5*(c'*z+c'*x + be'*lambda_E + bi'*lambda_I);
    condition = gap^2;
    condition = max(condition,max(tmp.*tmp));
    condition = max(condition,max(de.*de));
    condition = max(condition,max(xminusy.*xminusy));
    condition = max(condition,max(zminusy.*zminusy));
    
    iter = iter +1;
end
toc
