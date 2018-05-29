x = randn(10^8,1);

c = lines(5);

% Normal
figure; 
subplot(1,2,1); hold on;
binedges = -4:0.05:4;
Nbincenters = center(binedges);
[N,~] = histcounts(x,binedges); % Calculate histogram

% bar(Nbincenters,N./trapz(Nbincenters,N),1,'FaceColor',c(1,:))
% plot(Nbincenters,N./trapz(Nbincenters,N),'Color',c(1,:))
patch(Nbincenters,N./trapz(Nbincenters,N),c(1,:))
ylim([0,0.5])
plot([0,0],ylim,'Color',c(2,:))
plot([1,1].*mean(abs(x)),ylim,'Color',c(3,:))
ylabel('Probability Density')
xlabel('Sigma')
title('Normal distribution')

% Half-normal
subplot(1,2,2); hold on;
binedges = 0:0.05:4;
HNbincenters = center(binedges);
[HN,~] = histcounts(abs(x),binedges); % Calculate histogram

% bar(HNbincenters,HN./trapz(HNbincenters,HN),1,'FaceColor',c(1,:))
% plot(HNbincenters,HN./trapz(HNbincenters,HN),'Color',c(1,:))
patch([HNbincenters(1),HNbincenters],[0,HN./trapz(HNbincenters,HN)],c(1,:))
ylim([0,1])
plot([1,1].*mean(abs(x)),ylim,'Color',c(2,:))
ylabel('Probability Density')
xlabel('Sigma')
title('Half-Normal distribution')