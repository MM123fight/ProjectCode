function [b] = readVec(path,n)
    file = fopen(path);
    b_data = textscan(file,'%f');
    b = b_data{1};
    if(size(b) ~= n)
        error("vector size is not accordant in "+path);
    end
    fclose(file);
end