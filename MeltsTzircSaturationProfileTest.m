% cd ~/Desktop/meltstzirc/output/
cd ~/Desktop/zircon' cryst timescales'/

if ~exist('dist','var'); load zircondistribution; end;

%% draw from new PDF

nzircs=1;
agemax = 27.3;
agemin = 27.0;


agerange = agemax - agemin;


ages = agerange.*(1-randsample(xi,nzircs,true,dist)) + agemin;
ages = sort(ages)';
uncert = random('Gamma',5,0.02,size(ages))/1 .* ages;
observedages = ages + randn(size(ages)).*uncert;


sorted = sortrows([observedages, uncert],1);
observedages = sorted(:,1);
observeduncert = sorted(:,2);


figure; hold on;
plot([0,nzircs+1],[agemin agemin],'b');
plot([0,nzircs+1],[agemax agemax],'b');
errorbar(observedages,2*observeduncert,'.r')

%% Export and run image

exportmatrix(sorted,'zircondata.tsv','\t');

system('./tzircrystimage 200 MeltsTZircDistribtuion.csv zircondata.tsv > imagedata.tsv');


load imagedata.tsv

im = imagedata(2:end, 2:end);
x = imagedata(1,2:end);
y = imagedata(2:end,1);

figure; imagesc(x,y,im)

hold on; plot(agemax,agemin,'+k','MarkerSize',10,'LineWidth',1.5)
hold on; plot(max(observedages),min(observedages),'.k','MarkerSize',10)
hold on; plot(wmean(observedages,observeduncert),wmean(observedages,observeduncert),'.m','MarkerSize',10)
legend('True answer','First and last measured','Weighted  mean')
xlabel('Beginning of zircon crystallization')
ylabel('End of zircon crystallization')
title(['N = ' num2str(nzircs)]);
formatfigure 


%% Export and run metropolis sampler

exportmatrix(sorted,'zircondata.tsv','\t');
tic;
system('./tzircrystmetropolis 100 100000 2.9 MeltsTZircDistribtuion.csv zircondata.tsv > metropolisdata.tsv');
toc
load metropolisdata.tsv

figure; plot(metropolisdata(:,3))
figure; plot(metropolisdata(:,1:2))
transition_probability=sum(diff(metropolisdata(:,3))~=0)./length(metropolisdata(:,3))



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

%% Try and find original min and max
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