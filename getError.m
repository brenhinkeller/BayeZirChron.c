function [err_sigma, wmerr_sigma, yzerr_sigma] = getError(xi, dist, nzircs, dt_sigma, i)
agemax = 27.03;
agemin = 27.0;
agerange = agemax - agemin;

rng(sum(clock)*i);
ages = agerange.*(1-randsample(xi,nzircs,true,dist)) + agemin;
ages = sort(ages)';
uncert = agerange./dt_sigma.*ones(size(ages));
observedages = ages + randn(size(ages)).*uncert;

sorted = sortrows([observedages, uncert],1);
observedages = sorted(:,1);
observeduncert = sorted(:,2);

% Export and run metropolis sampler
workdir=['out' num2str(i)];
system(['mkdir -p ' workdir]);
exportmatrix(sorted,[workdir '/zircondata.tsv'],'\t');
system(['cd ' workdir '; ../tzircrystmetropolis 100000 ../MeltsTZircDistribtuion.tsv ./zircondata.tsv > metropolisdata.tsv']);
load([workdir '/metropolisdata.tsv']);

[wm, ~, ~]=wmean(observedages,observeduncert);

err_sigma = (nanmean(metropolisdata(10000:end,1))-agemin)/nanmean(observeduncert);
wmerr_sigma = (wm-agemin)/nanmean(observeduncert);
yzerr_sigma = (min(observedages)-agemin)/nanmean(observeduncert);
end
