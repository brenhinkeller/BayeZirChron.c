% cd ~/Desktop/meltstzirc/output/
cd ~/Desktop/zircon' cryst timescales'/

if ~exist('dist','var'); load zircondistribution; end;

%% draw from new PDF

nzircs=100;
agemax = 27.3;
agemin = 27.0;


agerange = agemax - agemin;


ages = agerange.*(1-randsample(xi,nzircs,true,dist)) + agemin;
ages = sort(ages)';
uncert = random('Gamma',5,0.02,size(ages))/100 .* ages;
observedages = ages + randn(size(ages)).*uncert;


sorted = sortrows([observedages, uncert],1);
observedages = sorted(:,1);
observeduncert = sorted(:,2);


figure; hold on;
plot([0,nzircs+1],[agemin agemin],'b');
plot([0,nzircs+1],[agemax agemax],'b');
% errorbar(ages,2*uncert,'.b')
errorbar(observedages,2*observeduncert,'.r')

%% Export and run

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