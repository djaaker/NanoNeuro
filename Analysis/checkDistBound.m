function result = checkArrayBound(array, bound)
    % Initialize result as true
    result = true;
    minValue = min(array);
    maxValue = max(array);

    if(maxValue-minValue<=bound)
        result = true;
    else
        result = false;
    end
    
end
