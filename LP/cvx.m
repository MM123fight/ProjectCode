root = "/Users/lumeng/Desktop/CCode/LP/data";
experiment_names = {"Mat100|100_1000_0.01_0.5","test3","Mat1000|0_10000_0.01_0.5","Mat0|1000_10000_0.01_0.5"};
experiment = root + "/" + experiment_names{1};
%svm_file = experiment + "/LPtoSVM";
%[b,A] = libsvmread(char(svm_file));
[nb,nf,mi,me,m,n,Ae,Ai,be,bi,c] = dataRead(experiment);



%{
[b,A] = libsvmread(char(svm_file+"/svmtoLP"));
[m,n] = size(A);
[c] = readVec(svm_file+"/c",n);
%}

bound = [zeros(nb,1);-Inf(nf,1)]
linprog(c',Ai,bi,Ae,be,bound); % optimize
c'*ans 


cvx_begin
   variable x(n)
   dual variables y z
   minimize( c' * x )
   subject to
      y: Ai * x <= bi
      z: Ae * x == be
      x(1:nb) >= 0
cvx_end  
cvx_optval 

Ab = [Ai(:,1:nb);Ae(:,1:nb)];
Af = [Ai(:,nb+1:n);Ae(:,nb+1:n)];
cb = c(1:nb);
cf = c(nb+1:n);
b = [bi;be];
cvx_begin
   variable x(n) 
   variable lambda(m)
   OBJ = 0.5*sum(square_pos(Ai * x- bi)) + 0.5*sum(square(Ae * x- be));
   OBJ = OBJ + 0.5*sum(square_pos(-Ab' * lambda - cb));
   OBJ = OBJ + 0.5*sum(square(-Af'*lambda -cf));
   OBJ = OBJ + 0.5*sum(square(c'*x + b'*lambda));
   minimize(OBJ)
   subject to
      x(1:nb) >= 0
      lambda(1:mi) >= 0
cvx_end
cvx_optval

%{
[m,n] = size(A);
lambda = 1;
cvx_begin
    variable x(n)
    minimize( lambda* sum(abs(x)) + 0.5*sum(square(A*x-b)))
cvx_end
cvx_optval
%}

%{
[m,n] = size(A);
mplus = mi;
nplus = nb+mi;
lambda = 0.;
if((mplus < m)&&(mplus >0))
    cvx_begin
    variable x(n)
    OBJ = 0.5*sum(square_pos(A(1:mplus,:)* x- b(1:mplus)));
    OBJ = OBJ + 0.5*sum(square(A(mplus+1:m,:)* x- b(mplus+1:m))) ;
    OBJ = OBJ + 0.5*lambda * sum(square(x));
    minimize(OBJ)
    subject to
    x(1:nplus) >= 0
    cvx_end
    cvx_optval
elseif(mplus == 0)
    cvx_begin
    variable x(n)
    minimize(0.5 * lambda* sum(square(x)) + 0.5*sum(square(A*x-b)))
    subject to
    x(1:nplus) >= 0
    cvx_end
    cvx_optval 
elseif(mplus == m)
    cvx_begin
    variable x(n)
    OBJ = 0.5 * lambda* sum(square(x)) + 0.5 * sum(square_pos(A(1:m,:)* x- b(1:m)));
    minimize(OBJ)
    subject to
    x(1:nplus) >= 0
    cvx_end
    cvx_optval   
end

%}
%{
[b_data,A_data] = libsvmread(char(experiment));
[N,d] = size(A_data);

c = ones(N,1);
m = (k-1)*N;
n1 = k*d;
n = n1+N;
value;
for i = 1:N
    for j = 1:kf
        if(j<d(i))
            value = [A(j,:);A(d(i),:)];
            
            
    end
end
%}