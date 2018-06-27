%% Plot MSWD as a function of N for different dt_sigma

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
    for j = 1:length(nr)
        mswdr(j) = nanmean(abs(mswd(n==nr(j))));
        mswdsigma(j) = nanstd(abs(mswd(n==nr(j))));
    end
    errorbar(nr, mswdr, 2*mswdsigma,'-');
end
set(gca,'Xscale','log')
xlabel('Number of zircons')
ylabel('MSWD')
xlim([min(nr) max(nr)]);
legend(cellfun(@(x) num2str(x,'Dt = %gsigma'), num2cell(dts),'UniformOutput',false))
xlim([9 1050]);
set(gca,'YTick',[0 1 2 3 4 5]);
% formatfigure;