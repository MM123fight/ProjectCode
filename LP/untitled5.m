%read data
root = "./data";
experiment_names = {"test3","Mat10|100_1000_0.1_0.5","Mat100|1000_10000_0.01_0.5"};
experiment = root + "/" + experiment_names{2};
[nb,nf,mi,me,m,n,Ae,Ai,be,bi,c] = dataRead(experiment);
%parameter
tol = 1e-6;
rho = 1;

%tau = 0.5*(1 + sqrt(5));
tau = 1;
sigma = tau*rho;
Ai = [Ai;Ae;-Ae];
bi = [bi;be;-be];
mi = 2*me+mi;

%primal variables
y = zeros(n,1);
z = zeros(n,1);
Xi = zeros(mi,1);
%dual variables
lambda_I = zeros(mi,1);
lambda_z = zeros(n,1);


tic;
t1 = clock;
AiT = Ai';
Invi =  inv(Ai*AiT+ eye(mi));
t2 = clock;
time = etime(t2,t1);

%ADMM
condition = 1;
iter = 1;
while(condition > tol)
    %{
    if(mod(iter,100) ==0)
        rho = rho+1;
    end
    %}
    
    t1 = clock;
    wi = lambda_I/rho -bi + Xi;
    vi = (lambda_z+c)/rho - y;
    z =  -vi +AiT*(Invi*(Ai*vi-wi));
 
    %step2:update y 
    y = z+ lambda_z/rho;
    y(1:nb) = max(y(1:nb),0);
    di = Ai*z - bi;
    tmp_vec =  di + lambda_I/rho;
    Xi = max(-tmp_vec,0);
    
    %step3:update lambda
    zminusy = z-y;
    tmp = di+Xi;
    lambda_I = lambda_I + rho*tmp;
    lambda_z = lambda_z + rho*zminusy;
    t2 = clock;
    time = time + etime(t2,t1);
     
    %check the stopping condition 
    optimal = c'*z;
    gap = optimal + (lambda_I'*bi);
    relative_gap = sqrt(gap^2/(optimal^2));
    condition = max(max(tmp.*tmp),max(zminusy.*zminusy));
    condition = max(condition, relative_gap);
    iter = iter +1;
    
end
toc