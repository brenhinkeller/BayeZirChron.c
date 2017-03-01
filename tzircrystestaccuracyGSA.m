v = 46;
cd(['~/Desktop/zircon cryst timescales/v' num2str(v) ' results/']);

% dts = [0.01 0.1 0.5 1 2 5 10 20 50 100];
dts = [0.01 0.5 1 2 5 10];

% dts = [0.01];

%% Mean absolute and relative error

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
    
    for j = 1:length(nr);
        mswdr(j) = nanmean(abs(mswd(n==nr(j))));
        datar(j,:) = nanmean(abs(data(n==nr(j),3:end)),1);
    end
    
    figure('OuterPosition',[300,300,1000,500]);
    subplot(1,2,1); hold on; plot(repmat(nr,1,5),datar(:,1:5));
    set(gca,'Xscale','log')
%     set(gca,'Yscale','log')
    xlabel('Number of zircons')
    ylabel('Average error (sigma)')
    title(['Delta t  =  ' num2str(dts(i)) ' sigma, MSWD: ' num2str(mswdr(end),'%.2f')])
    xlim([min(nr) max(nr)]);
    ylim([0 3])
    formatfigure
    legend({'Weighted Mean','MSWD Test','Youngest Zircon','Bayesian estimate','Bayesian (uniform prior)'},'Location','northwest')

        
    subplot(1,2,2); hold on; plot(repmat(nr,1,5), datar(:,6:10));
    plot(nr,ones(size(nr)),'--k');
    set(gca,'Xscale','log')
%     set(gca,'Yscale','log')
    xlabel('Number of zircons')
    ylabel('Average error / estimated error')
    title(['Delta t  =  ' num2str(dts(i)) ' sigma, MSWD: ' num2str(mswdr(end),'%.2f')])
    xlim([min(nr) max(nr)]);
    %     ylim([0.1 100])
    ylim([0 6])
    formatfigure
    
end


%%
v = 46;
cd(['~/Desktop/zircon cryst timescales/v' num2str(v) ' results/']);
% dts = [0.01 0.1 0.5 1 2 5 10 20 50 100];
dts = [5 2 0.5 0.01];

figure; hold on;
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
    mswdsigma = zeros(size(nr));
    for j = 1:length(nr);
        mswdr(j) = nanmean(abs(mswd(n==nr(j))));
        mswdsigma(j) = nanstd(abs(mswd(n==nr(j))));
    end
    errorbar(nr, mswdr, 2*mswdsigma,'-');
end
set(gca,'Xscale','log')
%     set(gca,'Yscale','log')
xlabel('Number of zircons')
ylabel('MSWD')
xlim([min(nr) max(nr)]);
legend(cellfun(@(x) num2str(x,'Dt = %gsigma'), num2cell(dts),'UniformOutput',false))
formatfigure;
xlim([9 1050]);
set(gca,'YTick',[0 1 2 3 4 5]);


%% Plot example datasets at different dt/sigma
cd ~/Desktop/zircon' cryst timescales'/
if ~exist('dist','var'); load zircondistribution; end;
dts = [0.01 0.5 1 2 5 10];
dts = 7;

nzircs = 10;
agemax = 100.1;
agemin = 100;
agerange = agemax - agemin;

for i=1:length(dts)
    dt_sigma =dts(i);
    
    ages = agerange.*(1-randsample(xi,nzircs,true,dist)) + agemin;
    ages = sort(ages)';
    % uncert = random('Gamma',5,0.02,size(ages))/100 .* ages;
    uncert = agerange./dt_sigma.*ones(size(ages));
    observedages = ages + randn(size(ages)).*uncert;
    
    
    sorted = sortrows([observedages, uncert],1);
    observedages = sorted(:,1);
    observeduncert = sorted(:,2);
    
    
    figure; hold on;
    errorbar(observedages,2*observeduncert,'.')
    plot([0,nzircs+1],[agemin agemin]);
    plot([0,nzircs+1],[agemax agemax]);
%    plot([0,nzircs+1],[agemin agemin],'b');
%    plot([0,nzircs+1],[agemax agemax],'b');
%    errorbar(observedages,2*observeduncert,'.r')
    xlabel('');
    xlim([0,nzircs+1]);
    set(gca,'Xtick',[]);
    title(['Delta t  =  ' num2str(dts(i)) ' sigma'])
    formatfigure
end


%% Plot example datasets at different N
cd ~/Desktop/zircon' cryst timescales'/
if ~exist('dist','var'); load zircondistribution; end;
Ns = [10 100];
dt_sigma = 0.01;

agemax = 100.1;
agemin = 100;
agerange = agemax - agemin;

for i=1:length(Ns)
    nzircs = Ns(i);
    
    ages = agerange.*(1-randsample(xi,nzircs,true,dist)) + agemin;
    ages = sort(ages)';
    % uncert = random('Gamma',5,0.02,size(ages))/100 .* ages;
    uncert = agerange./dt_sigma.*ones(size(ages));
    observedages = ages + randn(size(ages)).*uncert;
    
    
    sorted = sortrows([observedages, uncert],1);
    observedages = sorted(:,1);
    observeduncert = sorted(:,2);
    
    
    figure; hold on;
    errorbar(observedages,2*observeduncert,'.')
    plot([0,nzircs+1],[agemin agemin]);
    plot([0,nzircs+1],[agemax agemax]);
%    plot([0,nzircs+1],[agemin agemin],'b');
%    plot([0,nzircs+1],[agemax agemax],'b');
%    errorbar(observedages,2*observeduncert,'.r')
    xlabel('');
    xlim([0,nzircs+1]);
    set(gca,'Xtick',[]);
    title(['Delta t  =  ' num2str(dt_sigma) ' sigma, N = ' num2str(nzircs)])
    formatfigure
end