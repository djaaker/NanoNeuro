function result = checkArrayBound(array, lowerBound, upperBound)
    % Initialize result as true
    result = true;
    
    % Iterate over the array
    for i = 1:length(array)
        % Check if the current element is outside the bound
        if array(i) < lowerBound || array(i) > upperBound
            % Update result to false
            result = false;
            % Exit the loop
            break;
        end
    end
end
