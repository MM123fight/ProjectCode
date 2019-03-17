%read data matlab中文乱码，请先在记事本里输入中文，然后复制到matlab里
root = "./data";
experiment_names = {"colon","Mat0|100_1000_0.1_0.5","Mat0|1000_10000_0.01_0.5"};
experiment = root + "/" + experiment_names{1};
svm_file = experiment + "/SVMtoLP";
[nb,nf,mi,me,m,n,Ae,Ai,be,bi,c] = dataRead(svm_file);
%[nb,nf,mi,me,m,n,Ae,Ai,be,bi,c] = dataRead(experiment);


%parameter for outer iterations
rho = 1;
tol = 1e-3;
tau = 0.5*(1 + sqrt(5));
sigma = tau*rho;

%subtype = "NU_ACDM";
subtype = "exact";
ACDM_tol = tol;
cond_type = "grad_cord";
scale = 1;

%ADMM1
n_bar = n + mi;
Ae = [Ae zeros(me,mi);Ai eye(mi)];
c = [c;zeros(mi,1)];
be = [be;bi];
me = me + mi;

%primal variables
x = zeros(n_bar,1);
y = zeros(n_bar,1);
%dual variables
lambda_E = zeros(me,1);
lambda_x = zeros(n_bar,1);

tic;
t1 = clock;
AeT = Ae';
Inv =  inv(Ae*AeT+ eye(me));
t2 = clock;
time = etime(t2,t1);
%ADMM
condition = 1;
iter = 1;

while(condition > tol)
    t1 = clock;
    we = lambda_E/rho - be;
    ve = (lambda_x+c)/rho - y;
    if(subtype == "exact")
        Aeveminuswe = Ae*ve-we;
        g = Inv*Aeveminuswe;
        x = -ve +AeT*g;
    elseif(subtype == "NU_ACDM")

        ge = -(AeT*we+ve);
        x = NU_ACDM(n_bar,Ae,ge,scale,zeros(n_bar,1),ACDM_tol,cond_type);
    end

    %xsquare = x'*x
   
    %step2:update y
    y = x+lambda_x /rho;
    y(1:nb) = max(y(1:nb),0);
    y(n+1:n_bar) = max(y(n+1:n_bar),0);
    %ysquare = y'*y
    %step3:update lambda
    de = Ae*x - be;
    xminusy = x-y;
    lambda_E = lambda_E + rho*de;
    lambda_x = lambda_x + rho*xminusy;
    t2 = clock;
    time = time + etime(t2,t1);
     
    %check the stopping condition
    optimal = c'*x
    gap = optimal + (lambda_E'*be);
    relative_gap = sqrt(gap^2/(optimal^2));
    condition = max(max(de.*de),max(xminusy.*xminusy));
    condition = max(condition,relative_gap);
    iter = iter+1;
end
toc

