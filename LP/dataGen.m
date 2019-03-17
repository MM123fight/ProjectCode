
me = 100;
mi = 100;
nb = 990;
nf = 10;
n = nf + nb;
density = 0.01;
condition = 0.5;
experiment = "Mat" + num2str(me) + "|" + num2str(mi) + "_" + num2str(n)+"_" +num2str(density) +"_" + num2str(condition);

experiment_file = "/Users/lumeng/Desktop/CCode/LP/data/" + experiment;
mkdir(experiment_file);
mkdir(experiment_file+"/result");
meta_path = experiment_file + "/meta";
Ai_path = experiment_file + "/A";
Ae_path = experiment_file + "/Aeq";
bi_path = experiment_file + "/b";
be_path = experiment_file + "/beq";
c_path = experiment_file + "/c";
svm_path = experiment_file + "/svm";

Ai = sparse(mi,n);
bi = zeros(mi,1);
Ae = sparse(me,n);
be = zeros(me,1);
x = ones(n,1);
c = rand(n,1);
if(mi > 0)
    Ai = sprand(mi,n,density,condition);
    bi = Ai*x +1;
end

if(me > 0)
    Ae = sprand(me,n,density,condition);
    be = Ae*x;
end

writemeta(meta_path,nb,nf,mi,me);
writeMat(Ai_path,Ai,mi,n);
writeMat(Ae_path,Ae,me,n);
writeVec(bi_path,bi);
writeVec(be_path,be);
writeVec(c_path,c);
b = [bi;be];
A =[Ai;Ae];
libsvmwrite(char(svm_path),b, A);


function  writemeta(path,nb,nf,mi,me)
    name = ["nb","nf","mi","me"];
    size = [nb,nf,mi,me];
    file = fopen(path,'w');
    for row = 1:3
        fprintf(file, '%s %d\n', name(row), size(row));
    end 
    fclose(file);
    row = 4;
    file = fopen(path,'a');
    fprintf(file, '%s %d', name(row), size(row));
    fclose(file);
end
function writeMat(path,A,m,n)
    file = fopen(path,'w');
    if(m == 0)
        fprintf(file,'%d %d %.1f',m,n,0);
        fclose(file);
    else
        fprintf(file,'%d %d %.1f\n',m,n,0);
        [col,row,val] = find(A');
        s = length(row);
        if(s==1)
            fprintf(file,'%d %d %f',[row(s);col(s);val(s)]);
            fclose(file);
        else
            fprintf(file,'%d %d %f\n',[row(1:s-1)';col(1:s-1)';val(1:s-1)']);
            fclose(file);
            file = fopen(path,'a');
            fprintf(file,'%d %d %f',[row(s);col(s);val(s)]);
            fclose(file);
        end
        
    end
    
  
end
function writeVec(path,b)
    file = fopen(path,'w');
    s = length(b);
    if(s == 1)
        fprintf(file, '%f', b);
        fclose(file);
    elseif(s>1)
        fprintf(file, '%f\n', b(1:s-1));
        fclose(file);
        file = fopen(path,'a');
        fprintf(file, '%f', b(s));
        fclose(file);
    end 
end

