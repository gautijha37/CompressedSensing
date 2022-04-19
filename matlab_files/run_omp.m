function [y] = run_omp(t, theta, observed_sig)
    maxsnr = 0;
    min_hr = 4/6;
    max_hr = 12/6;
    sparsity_range = 2 * t(end)/min_hr : 5 * t(end)/max_hr;
    for i = sparsity_range
        y1 =  algo_omp(i, theta, observed_sig);
end

