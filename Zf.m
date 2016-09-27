function out = Zf(f, mswd)
x = f./2;
out = exp(x.*log(x) - stirling_lgamma(x) + (x-1).*log(mswd) - x.*mswd); % Distribution of the MSWD from Wendt and Carl (1991)


