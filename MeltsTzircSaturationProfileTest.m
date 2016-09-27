% cd ~/Desktop/meltstzirc/output/
cd ~/Desktop/zircon' cryst timescales'/

if ~exist('dist','var'); load volcaniczircondistribution.mat; end;

%% draw from new PDF

nzircs = 1000;
dt_sigma =0.1;

agemax = 100.1;
agemin = 100;
agerange = agemax - agemin;

ages = agerange.*(1-randsample(xi,nzircs,true,dist)) + agemin;
ages = sort(ages)';
% uncert = random('Gamma',5,0.02,size(ages))/100 .* ages;
uncert = agerange./dt_sigma.*ones(size(ages));
observedages = ages + randn(size(ages)).*uncert;


sorted = sortrows([observedages, uncert],1);
observedages = sorted(:,1);
observeduncert = sorted(:,2);


figure; hold on;
plot([0,nzircs+1],[agemin agemin],'b');
plot([0,nzircs+1],[agemax agemax],'b');
errorbar(observedages,2*observeduncert,'.r')

[wm, wsigma, mswd]=wmean(observedages,observeduncert)
%Zf(nzircs-1,mswd)
Pf = Zf(nzircs-1,mswd)./Zf(nzircs-1,1)


% Export and run metropolis sampler

exportmatrix(sorted,'zircondata.tsv','\t');
system('./tzircrystmetropolis 1000 100000 2.9 MeltsTZircDistribtuion.csv zircondata.tsv > metropolisdata.tsv');
load metropolisdata.tsv

% figure; plot(metropolisdata(:,3))
figure; plot(metropolisdata(:,1:2))
transition_probability=sum(diff(metropolisdata(:,3))~=0)./length(metropolisdata(:,3))
averages = nanmean(metropolisdata(10000:end,1:2))
yz = min(observedages)
err_sigma = (nanmean(metropolisdata(10000:end,1))-agemin)/nanmean(observeduncert)
err_wmerr = (nanmean(metropolisdata(10000:end,1))-agemin)/(wm-agemin)
err_yzerr = (nanmean(metropolisdata(10000:end,1))-agemin)/(min(observedages)-agemin)

%% Export and run image

exportmatrix(sorted,'zircondata.tsv','\t');

system('./tzircrystimage 200 MeltsTZircDistribtuion.csv zircondata.tsv > imagedata.tsv');


load imagedata.tsv

im = imagedata(2:end, 2:end);
x = imagedata(1,2:end);
y = imagedata(2:end,1);

figure; imagesc(x,y,im)

hold on; plot(agemax,agemin,'+k','MarkerSize',10,'LineWidth',1.5)
% hold on; plot(max(observedages),min(observedages),'.k','MarkerSize',10)
hold on; plot(wmean(observedages,observeduncert),wmean(observedages,observeduncert),'+m','MarkerSize',10,'LineWidth',1.5)
legend('True answer','Weighted  mean')
xlabel('Beginning of zircon crystallization')
ylabel('End of zircon crystallization')
title(['N = ' num2str(nzircs)]);
formatfigure 

%% Check optimal step factor
stepfactor = (1.9:0.1:3)';
transition_probability=NaN(size(stepfactor));

pool=gcp; %Start a parellel processing pool if there isn't one already
parfor i=1:length(stepfactor);
    transition_probability(i) = getTransitionProbability(stepfactor(i))
%     system(['./tzircrystmetropolis 100 100000 ' num2str(stepfactor(i)) ' MeltsTZircDistribtuion.csv zircondata.tsv > metropolisdata.tsv']);
%     load metropolisdata.tsv
%     transition_probability(i) = sum(diff(metropolisdata(:,3))~=0)./length(metropolisdata(:,3));
end

figure; plot(stepfactor,transition_probability,'.','MarkerSize',20)

%% Check distribution of Zf for different N and dt/sigma
nsims = 1000;

nzircs = 100;
dt_sigma = 2;

agemax = 27.03;
agemin = 27.0;
agerange = agemax - agemin;
Zfs = NaN(nsims, 1);
mswds = NaN(nsims, 1);
for i = 1:nsims
    ages = (agerange.*(1-randsample(xi,nzircs,true,dist)) + agemin)';
    uncert = agerange./dt_sigma.*ones(size(ages));
    observedages = ages + randn(size(ages)).*uncert;
    [~, ~, mswds(i)]=wmean(observedages,uncert);
    Zfs(i) = Zf(nzircs-1,mswds(i));
end

mswds(mswds<1) = 1;

figure; hold on;
subplot(1,3,1); hist(Zfs,50); xlabel('Zf')
subplot(1,3,2); hist(Zfr(nzircs-1,mswds,1),50); xlabel('Zf/Zfmax')
% subplot(1,3,2); hist(Zfs./Zf(nzircs-1,1-2/(nzircs-1)),50); xlabel('Zf/Zfmax')
subplot(1,3,3); hist(mswds,50); xlabel('MSWD')
% subplot(1,3,3); hist(Zfs./Zf(nzircs-1,1)./mswds,50); xlabel('Zf/Zfmax/Mswd')
% subplot(1,3,3); hist(Zfr(nzircs-1,mswds,1)./mswds,50); xlabel('Zf/Zfmax/Mswd')


%% Explore average errors

nzircs = 7;
dt_sigma = 0.01;
nsims = 48;

err_sigma = NaN(nsims,1);
wmerr_sigma = NaN(nsims,1);
yzerr_sigma = NaN(nsims,1);
gcp;
parfor i=1:nsims
    [err_sigma(i), wmerr_sigma(i), yzerr_sigma(i)] = getError(xi, dist, nzircs, dt_sigma, i)
end

% err_sigma
% wmerr_sigma
% yzerr_sigmaer

fprintf('Averages\n');
err_sigma_ave = nanmean(err_sigma)
wmerr_sigma_ave = nanmean(wmerr_sigma)

yzerr_sigma_ave = nanmean(yzerr_sigma)

fprintf('Absolute Averages\n');
err_sigma_abs = nanmean(abs(err_sigma))
wmerr_sigma_abs = nanmean(abs(wmerr_sigma))
yzerr_sigma_abs = nanmean(abs(yzerr_sigma))

fprintf('\n');


%% Try and find original min and max (old matlab version)
% 
% nsims = 10000;
% 
% imax = NaN(nsims,1);
% imin = NaN(nsims,1);
% misfit = NaN(nsims,1);
% 
% imin(1) = min(observedages);
% imax(1) = max(observedages);
% 
% plot([0,nzircs+1],[imin(1) imin(1)],'r');
% plot([0,nzircs+1],[imax(1) imax(1)],'r');
% 
% 
% iages = sort(((imax(1)-imin(1)).*(1-randsample(xi,nzircs,true,dist)) + imin(1)))';
% misfit(1) = sum(((observedages - iages)./observeduncert).^2);
% 
% lastiminchange = imax(1)-imin(1);
% lastimaxchange = imax(1)-imin(1);
% 
% for i=2:nsims
%     if rand>0.4
%         imin(i) = imin(i-1) + randn.*2.*lastiminchange;
%         imax(i) = imax(i-1);
%         
%         iages = sort(((imax(i)-imin(i)).*(1-randsample(xi,nzircs,true,dist)) + imin(i)))';
%         misfit(i) = sum(((observedages - iages)./observeduncert).^2);
%         
%         
%         if misfit(i)<misfit(i-1)
%             % Yay!
%             lastiminchange = imin(i)-imin(i-1);
%         else
%             imin(i) = imin(i-1);
%             imax(i) = imax(i-1);
%             misfit(i) = misfit(i-1);
%         end
%     else
%         imin(i) = imin(i-1);
%         imax(i) = imax(i-1) + randn.*2.*lastimaxchange;
%         
%         iages = sort(((imax(i)-imin(i)).*(1-randsample(xi,nzircs,true,dist)) + imin(i)))';
%         misfit(i) = sum(((observedages - iages)./observeduncert).^2);
%         
%         
%         if misfit(i)<misfit(i-1)
%             % Yay!
%             lastimaxchange = imax(i)-imax(i-1);
%         else
%             imin(i) = imin(i-1);
%             imax(i) = imax(i-1);
%             misfit(i) = misfit(i-1);
%         end
%     end
% end
% 
% 
% plot([0,nzircs+1],[imin(i) imin(i)],'m');
% plot([0,nzircs+1],[imax(i) imax(i)],'m');
% 
% 
% figure; plot(imin,'r')
% hold on; plot(imax,'b')
% 
% 
%     