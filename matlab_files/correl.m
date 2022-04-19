function [correlation_vec] = correl(A)
[~, n] = size(A);
% correlation_vec = zeros(1, n);

[~,correlation_vec] = max(abs(A' * A));

end

