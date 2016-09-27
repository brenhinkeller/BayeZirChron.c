function out = Zfr(f, mswd, origin)
out = exp((f./2-1).*log(mswd) - f./2.*(mswd-origin));
% out(mswd<1)=1;