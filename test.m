nzircs = 3;
nsims = 100;
wxs = NaN(nsims,1);
wsigmas = NaN(nsims,1);
true = 100;
zsigma = 1;



zircsigmas = zsigma*ones(nzircs,1);
for i=1:nsims
    zircs = true + randn(nzircs,1)*zsigma;
    [wxs(i),wsigmas(i),~] = wmean(zircs, zircsigmas);
end


nanmean(abs(wxs-true))*1.253
nanmean(wsigmas)

nanmean(abs(wxs-true))*1.253/nanmean(wsigmas)