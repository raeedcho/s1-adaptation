function [mean_tensor,sigma] = extract_tme_features(data_tensor)
    % Function to extract statistical features from firing rate tensor, including:
    %   - full mean tensor
    %   - covariances for each condition (sigma{<tensor dim>})
    % For the curl field task, tensor dimensions are:
    %   neurons(N) x time (T) x target_block (C) x learning_block (L)
    % Based on extractFeatures.m in TME library
    tensor_shape = size(data_tensor);
    sigma = cell(1,length(tensor_shape));

    % calculate mean tensor (three possible methods)
    %% Mean is caluclated by subtracting each reshaping mean
    % generic reshaping mean subtraction
    data_tensor0 = data_tensor;
    tensor_dims = 1:length(tensor_shape);
    for dimnum = 1:length(tensor_dims)
        sum_dim = tensor_dims(tensor_dims~=dimnum);
        dim_mean = sumTensor(data_tensor0, sum_dim)/prod(tensor_shape(sum_dim));
        data_tensor0 = bsxfun(@minus, data_tensor0,dim_mean);
    end
    mean_tensor = data_tensor-data_tensor0;

    %% subtract the mean tensor and calculate the covariances
    X0 = data_tensor-mean_tensor;
    for dimnum = 1:length(tensor_dims)
        sum_dim = tensor_dims(tensor_dims~=dimnum);
        permutation = [sum_dim dimnum];
        X_reshaped = reshape(permute(X0,permutation),[],tensor_shape(dimnum));
        sigma{dimnum} = (X_reshaped'*X_reshaped);
    end
end

