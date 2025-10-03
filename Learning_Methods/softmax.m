function y = softmax(x)
    % Compute softmax function
    exp_x = exp(x - max(x));  % Subtract max for numerical stability
    y = exp_x / sum(exp_x);
end