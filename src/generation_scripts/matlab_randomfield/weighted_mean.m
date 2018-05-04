function [average] = weighted_mean(X,weights)
    if ~all(size(weights)==size(X))
        warning('Size of data does not match size of weights.');
    end
    X_vec = reshape(X,[1,size(X,1)*size(X,2)*size(X,3)]);
    weights_vec = reshape(weights,[1,size(weights,1)*size(weights,2)*size(weights,3)]);
    average = sum(weights_vec.*X_vec)./sum(weights_vec);
end

