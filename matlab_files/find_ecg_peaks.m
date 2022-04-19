function [sparsity] = find_ecg_peaks(sig)
mx = max(sig);
sparsity = 0;
for i = 1 : length(sig)
    if(sig(i) >= 0.9*mx)
        sparsity = sparsity + 1;
    end
end

