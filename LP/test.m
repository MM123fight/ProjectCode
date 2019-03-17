root = "./data";
experiment_names = {"test1","test2","Mat10|100_1000_0.1_0.5"};
experiment = root + "/" + experiment_names{3};
[nb,nf,mi,me,m,n,Ae,Ai,be,bi,c] = dataRead(experiment);
x = zeros(n,1);
w = Ai*x -bi;
v = x + c;
d = ones(n,1);
tol = 1e-6;
rho = 1;
%x = SSN(mi,n,Ai,x,w,v,tol,c,"exact");
%x1 = (Ai'*Ai+eye(n))\d;

%x = CG_inv(mi,Ai,d,tol,1);


%{
tic
for i=1:10
    x0 = d;
    x =CG(Ai,d,x0,tol,1);
end
toc
%}



x0 = zeros(n,1);
A = Ai;
x_star = ones(n,1);
d = ones(n,1);
scale = 1;
tol = 1e-3;
cond_type = "grad_cord";
tic
x = NU_ACDM(n,A,d,scale,x0,tol,cond_type)
toc
v = (A'*A +eye(n))*x - d;
v'*v
%}
%{

tol = 1e-6;
z0 = zeros(n,1);

rho = 1;
lambda_I = zeros(mi,1);
lambda_z = zeros(n,1);
y = zeros(n,1);

w = Ai*z0+lambda_I/rho-bi;
v = z0+(lambda_z +c)/rho-y;
 
cond_type = "grad_sum";
tic
 z = SSN(mi,n,Ai,z0,w,v,tol,c,rho);
toc
%}


