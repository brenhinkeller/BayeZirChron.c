%% Generate fake zircon data

nsamples=12;
eruptionage=rand * 100;
agerange=1.3;
ageuncert=0.1;
ageuncerthet=2;

trueages=rand(nsamples,1).*agerange+eruptionage;
lastzircon=min(trueages);
measureduncert=ageuncert .* (rand(nsamples,1).*ageuncerthet +1);
measuredages=trueages + randn(nsamples,1) .* measureduncert;

[measuredages, IX] = sort(measuredages);
measureduncert = measureduncert(IX);
trueages = trueages(IX);

figure; errorbar(measuredages, measureduncert*2,'.')
hold on; plot([0, nsamples],[lastzircon, lastzircon],'c')
hold on; plot([0, nsamples],[eruptionage, eruptionage],'r')

hold on; plot([0, nsamples],[eruptionage+agerange, eruptionage+agerange],'b')
