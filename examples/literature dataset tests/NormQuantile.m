function q = NormQuantile(F)
    % How far away from the mean (in units of sigma) should we expect proportion
    % F of the samples to fall in a Normal (Gaussian) distribution
    q = sqrt(2)*erfinv(2*F-1);
end
