% Note: Make sure we're in the right folder before running! should be in
% /examples/synthetic dataset tests/

% Range of dt/sigma to explore
dts = [0.01 1 2 10];


%% Plot mean absolute error for each dt_sigma

% figure('OuterPosition',[600,0,400,1024]);

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
    nr = unique(n);
    nr = nr(arrayfun(@(x) sum(n==x)>10,nr));
    
    mswdr = zeros(size(nr));
    datar = zeros(length(nr),10);
    
    for j = 1:length(nr)
        mswdr(j) = nanmean(abs(mswd(n==nr(j))));
        datar(j,:) = nanmean(abs(data(n==nr(j),3:end)),1);
    end
    
%     subplot(length(dts),1,i);
    figure;
    plot(repmat(nr,1,5),datar(:,1:5));
    set(gca,'Xscale','log')
    xlabel('Number of zircons')
    ylabel('Average error (sigma)')
    title(['Delta t  =  ' num2str(dts(i)) ' sigma, MSWD: ' num2str(mswdr(end),'%.2f')])
    xlim([min(nr) max(nr)]);
    ylim([0 3])
    legend({'Weighted Mean','MSWD Test','Youngest Zircon','Bayesian estimate','Bayesian (uniform prior)'},'Location','northwest')
    formatfigure;
    fig = gcf;
    fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
    saveas(fig,['dt=', num2str(dts(i)), 'sigmaAbs.pdf'])
end

%% Plot mean relative error for each dt_sigma
% figure('OuterPosition',[600,0,400,1024]);

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
    nr = unique(n);
    nr = nr(arrayfun(@(x) sum(n==x)>10,nr));
    
    mswdr = zeros(size(nr));
    datar = zeros(length(nr),10);
    
    for j = 1:length(nr)
        mswdr(j) = nanmean(abs(mswd(n==nr(j))));
        datar(j,:) = nanmean(abs(data(n==nr(j),3:end)),1);
    end
    
%     subplot(length(dts),1,i); hold on;
    figure; hold on;
    plot(repmat(nr,1,5), datar(:,6:10));
    plot(nr,ones(size(nr)),'--k');
    set(gca,'Xscale','log')
    xlabel('Number of zircons')
    ylabel('Average error / estimated error')
    title(['Delta t  =  ' num2str(dts(i)) ' sigma, MSWD: ' num2str(mswdr(end),'%.2f')])
    xlim([min(nr) max(nr)]);
    ylim([0 6])
    formatfigure;
    fig = gcf;
    fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
    saveas(fig,['dt=', num2str(dts(i)), 'sigmaRel.pdf'])
end
