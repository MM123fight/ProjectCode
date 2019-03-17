%read data
root = "./data";
experiment_names = {"test3","Mat0|100_1000_0.1_0.5"};
experiment = root + "/" + experiment_names{1};
[nb,nf,mi,me,m,n,Ae,Ai,be,bi,c] = dataRead(experiment);

%parameter
tol = 1e-15;
rho = 1;

%ADMM1
n_bar = n + mi;
A = [Ae zeros(me,mi);Ai eye(mi)];
b = [be;bi];
AT = A';


y = zeros(n_bar,1);
z = zeros(m+n_bar,1);
g = zeros(n_bar,1);
H = A*AT + eye(m);
condition = 0.1;
tic
while(condition > tol)
    g(1:n) = (c+z(m+1:n+m))/rho -y(1:n);
    g(n+1:n_bar) = z(m+n+1:m+n_bar)/rho - y(n+1:n_bar);
    h = A*g - z(1:m)/rho + b;
    x = AT*(H\h) - g;
    y = z(m+1:m+n_bar)/rho +x;
    y(1:nb) = max(y(1:nb),0);
    d = A*x - b;
    xminusy = x - y;
    z = z+rho*[d;xminusy];
    condition = max(max(d.*d),max(xminusy.*xminusy));
end
toc
