experiment = "/Users/lumeng/Desktop/CCode/LP/data/test3";
[nb,nf,mi,me,m,n,Ae,Ai,be,bi,c] = dataRead(experiment);
mplus = mi
tau = 1;
theta = tau/n;
A = [Ai;Ae];
full(A);
b = [bi;be];
d = zeros(n,1);
y = zeros(n,1);
x = zeros(n,1);
z = x;
alpha = 1;

eta = (tau-1)/max(1,n-1);
nnz = sum(A~=0,2);
beta = 1 + eta*(nnz-1);
beta_sum = sum(beta);
v = sum(beta.*(A.*A));
v = v + alpha;
%K = ceil(2 * exp(1) * n * (sqrt(1+ 1./mu) -1)/tau + 1);
cord = [2;1;1;2;1;1;2;2;2;2;2;2;1;2;2;1;1;1;];
%cord=[3;1;1;3;1;2;3;3;2;3;3;2;1;3;2;1;2;2;1;1;1;2;1;3;2;3;2;2;2;1;1;1;1;2];
%cord = [9998;1630;2827;9473;2317;4850;9575;7444];
for i = 1:length(cord)
    theta_square = theta * theta;
    s = theta_square * y + z;
    w_sum = A*s-b;
    w_sum(w_sum(1:mplus) <0) = 0;
    grad = A(:,cord(i))'*w_sum + alpha*(s(cord(i))-d(cord(i)))
    t = - (grad * tau)/(n*theta*v(cord(i)));
    if((z(cord(i))+ t < 0)&& (cord(i) <=nb))
        t = -z(cord(i));
        z(cord(i)) = 0;
    else
        z(cord(i)) = z(cord(i))+ t;
    end
    t
    y(cord(i)) = y(cord(i))- (1-n*theta/tau)/theta_square*t;
    theta = 0.5*(sqrt(theta_square*theta_square + 4*theta_square) - theta_square);
end
    




