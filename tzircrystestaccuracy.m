cd ~/Desktop/zircontimescales/

% dts = [0.1 0.5 1 2 5 10 20];
dts = [0.1 0.5 1 5 10 20];

figure; 
for i=1:length(dts)
    
    data = load(['eruptionestimates' num2str(dts(i)) '.tsv']);
    
    n = data(:,1);
    wm = data(:,2);
    mt = data(:,3);
    yz = data(:,4);
    bz = data(:,5);
    bm = data(:,6);
    
    nr = unique(n);
    wmr = zeros(size(nr));
    mtr = zeros(size(nr));
    yzr = zeros(size(nr));
    bzr = zeros(size(nr));
    bmr = zeros(size(nr));
    
    for j = 1:length(nr);
        wmr(j) = nanmean(wm(n==nr(j)));
        mtr(j) = nanmean(mt(n==nr(j)));
        yzr(j) = nanmean(yz(n==nr(j)));
        bzr(j) = nanmean(bz(n==nr(j)));
        bmr(j) = nanmean(bm(n==nr(j)));
  
    end
    
    subplot(3,2,i); hold on;
    plot(nr, wmr,'r','LineWidth',1);
    plot(nr, mtr,'m','LineWidth',1);
    plot(nr, yzr,'c','LineWidth',1);
    plot(nr, bzr,'b','LineWidth',1);
    plot(nr, bmr,'k','LineWidth',1);
    
    set(gca,'Xscale','log')
    xlabel('Number of zircons')
    ylabel('Average error (sigma)')
    title(['Delta t  =  ' num2str(dts(i)) ' sigma'])
    
end

legend('Weighted Mean','MSWD Test','Youngest Zircon','Metropolis','Metropolis+MSWD')
