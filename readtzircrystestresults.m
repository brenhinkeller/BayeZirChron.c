
nNs = [1 2 5 10 20];

% figure; subplot(2,2,1);
for j = 1:length(nNs);
    N = nNs(j);
    
    % Load data
    data = load(['eruptionestimates' num2str(N) '.tsv']);
    n = data(:,1);
    wm = data(:,2);
    mt = data(:,3);
    yz = data(:,4);
    bz = data(:,5);
    bm = data(:,6);
    
    % Find averages for each N
    nr = unique(n);
    wmr=zeros(size(nr));
    mtr=zeros(size(nr));
    yzr=zeros(size(nr));
    bzr=zeros(size(nr));
    bmr=zeros(size(nr));
    for i = 1:length(nr)
        wmr(i) = nanmean(wm(n==nr(i)));
        mtr(i) = nanmean(mt(n==nr(i)));
        yzr(i) = nanmean(yz(n==nr(i)));
        bzr(i) = nanmean(bz(n==nr(i)));
        bmr(i) = nanmean(bm(n==nr(i)));
    end
    
    % Plot averages
%     subplot(2,2,j); hold on;
    figure; hold on
    plot(nr,wmr,'r');
    plot(nr,mtr,'m');
    plot(nr,yzr,'c');
    plot(nr,bzr,'b');
    plot(nr,bmr,'k');
    set(gca,'XScale','log');
    
    xlabel('Number of zircons')
    ylabel('Average misfit / sigma')
    title(['dt = ' num2str(N) ' sigma'])
    formatfigure
    
    
end

legend('Weighted Mean','MSWD Test','Youngest Zircon','Bayesian Estimate')
