function result = check_me(me)
    %the size affordable for inverse
    if(me < 1000)
        result = true;
    else
        result = false;
    end
end