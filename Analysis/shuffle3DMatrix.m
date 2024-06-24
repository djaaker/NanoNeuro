function shuffled_matrix = shuffle3DMatrix(matrix3D, shuffle_dimension)
    % Get the size of the input matrix
    matrix_size = size(matrix3D);
    
    % Check if the specified dimension is valid
    if shuffle_dimension < 1 || shuffle_dimension > ndims(matrix3D)
        error('Invalid shuffle dimension');
    end
    
    % Create a permutation of indices for the specified dimension
    permuted_indices = randperm(matrix_size(shuffle_dimension));
    
    if shuffle_dimension == 3
        shuffled_matrix = matrix3D(:,:,permuted_indices);
    end
    if shuffle_dimension == 2
        shuffled_matrix = matrix3D(:,permuted_indices,:);
    end
    if shuffle_dimension == 1
        shuffled_matrix = matrix3D(permuted_indices,:,:);
    end

end