cd '~/Desktop/zircon cryst timescales/v5 results/'

% dts = [0.1 0.5 1 2 5 10 20];
% dts = [0.1 0.5 1 5 10 20];
dts = [0.01 0.1 0.5 1 2 5 10 50];
% dts = [0.01 0.5 2 5];
% dts = [50];
%%

figure; hold on
for i=1:length(dts)
    
    name = ['eruptionestimates' num2str(dts(i))];
    if ~exist([name '.tsv'],'file')
        % Parse log file to remove any line that doesn't have 12
        % tab-delimited entries (numerical, nan, or inf).
        system(['grep -e ''^[-0-9\.naif][0-9\.+naief]*\(\t[-0-9\.naif][0-9\.+naief]*\)\{11\}$'' ' name '.log > ' name '.tsv']);
    end
    data = load([name '.tsv']);
    
    n = data(:,1);
    mswd = data(:,2);
    wm = data(:,3);
    mt = data(:,4);
    yz = data(:,5);
    bz = data(:,6);
    bm = data(:,7);
    
    nr = unique(n);
    mswdr = zeros(size(nr));
    wmr = zeros(size(nr));
    mtr = zeros(size(nr));
    yzr = zeros(size(nr));
    bzr = zeros(size(nr));
    bmr = zeros(size(nr));
    
    for j = 1:length(nr);
        mswdr(j) = nanmean(mswd(n==nr(j)));
        wmr(j) = nanmean(wm(n==nr(j)));
        mtr(j) = nanmean(mt(n==nr(j)));
        yzr(j) = nanmean(yz(n==nr(j)));
        bzr(j) = nanmean(bz(n==nr(j)));
        bmr(j) = nanmean(bm(n==nr(j)));
  
    end
    
    subplot(4,2,i); hold on;
    plot(nr, mswdr,'LineWidth',1);
    plot(nr, wmr,'r','LineWidth',1);
    plot(nr, mtr,'m','LineWidth',1);
    plot(nr, yzr,'c','LineWidth',1);
    plot(nr, bzr,'b','LineWidth',1);
    plot(nr, bmr,'k','LineWidth',1);
    
    set(gca,'Xscale','log')
    xlabel('Number of zircons')
    ylabel('Average error (sigma)')
    title(['Delta t  =  ' num2str(dts(i)) ' sigma'])
    xlim([min(nr) max(nr)]);
end
legend('MSWD','Weighted Mean','MSWD Test','Youngest Zircon','Metropolis','Metropolis+MSWD')
% legend('Weighted Mean','MSWD Test','Youngest Zircon','Metropolis','Metropolis+MSWD')
%%
figure; hold on
for i=1:length(dts)
    
    name = ['eruptionestimates' num2str(dts(i))];
    if ~exist([name '.tsv'],'file')
        % Parse log file to remove any line that doesn't have 12
        % tab-delimited entries (numerical, nan, or inf).
        system(['grep -e ''^[-0-9\.naif][0-9\.\+naief]*\(\t[-0-9\.naif][0-9\.\+naief]*\)\{11\}$'' ' name '.log > ' name '.tsv']);
    end
    data = load([name '.tsv']);
        
    n = data(:,1);
    mswd = data(:,2);
    wm = data(:,8);
    mt = data(:,9);
    yz = data(:,10);
    bz = data(:,11);
    bm = data(:,12);
    
    nr = unique(n);
    mswdr = zeros(size(nr));
    wmr = zeros(size(nr));
    mtr = zeros(size(nr));
    yzr = zeros(size(nr));
    bzr = zeros(size(nr));
    bmr = zeros(size(nr));
    
    for j = 1:length(nr);
        mswdr(j) = nanmean(mswd(n==nr(j)));
        wmr(j) = nanmean(wm(n==nr(j)));
        mtr(j) = nanmean(mt(n==nr(j)));
        yzr(j) = nanmean(yz(n==nr(j)));
        bzr(j) = nanmean(bz(n==nr(j)));
        bmr(j) = nanmean(bm(n==nr(j)));
  
    end
    
    subplot(4,2,i); hold on;
%     plot(nr, mswdr,'LineWidth',1);
    plot(nr, wmr,'r','LineWidth',1);
    plot(nr, mtr,'m','LineWidth',1);
    plot(nr, yzr,'c','LineWidth',1);
    plot(nr, bzr,'b','LineWidth',1);
    plot(nr, bmr,'k','LineWidth',1);
    
    set(gca,'Xscale','log')
    xlabel('Number of zircons')
    ylabel('Average error / estimated error')
    title(['Delta t  =  ' num2str(dts(i)) ' sigma'])
    xlim([min(nr) max(nr)]);
end

% legend('MSWD','Weighted Mean','MSWD Test','Youngest Zircon','Metropolis','Metropolis+MSWD')
legend('Weighted Mean','MSWD Test','Youngest Zircon','Metropolis','Metropolis+MSWD')
