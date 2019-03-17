function result = check_J(Jsize,mi)
    if(Jsize < max(sqrt(mi),1000))
        result = true;
    else
        result = false;
    end
end
    