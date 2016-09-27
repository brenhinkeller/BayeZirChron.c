function out = Pf(f,mswd)
datarows = f+1;
if mod(datarows,2)
    out = 0;
    for v=0:round((datarows-1)./2-1)
        out = out + (f./2).^v ./ factorial(v) .* mswd.^v;
    end
    out = out .* exp(-f.*mswd./2);
else
    out = 0;
    for v=1:round(datarows./2-1)
        out = out + f.^(v-1) .* 2.^v .* factorial(v)/factorial(2.*v) .* mswd.^(v-0.5);
    end
    out = 1 - 2.*(sqrt(erf(f.*mswd)) - sqrt(f./2.*pi).*exp(-f.*mswd./2).*out);
end