function out = stirling_lgamma(x)
% Stirling's approximation for the gamma function
out = (x - 0.5)*log(x) - x + 0.5*log(2*pi) + 1/(12*x) + 1/(360 * x * x * x);
out(x<1) = 0;
