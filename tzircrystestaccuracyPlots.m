cd('~/Desktop/Bayesian cryst timescales/results/');

% dts = [0.01 0.1 0.5 1 2 5 10 20 50 100];
dts = [0.01 0.5 1 5];


%% Plot mean absolute error for each dt_sigma

% figure;
figure('OuterPosition',[600,0,400,1024]);

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
    
    subplot(length(dts),1,i); plot(repmat(nr,1,5),datar(:,1:5));
    set(gca,'Xscale','log')
    xlabel('Number of zircons')
    ylabel('Average error (sigma)')
    title(['Delta t  =  ' num2str(dts(i)) ' sigma, MSWD: ' num2str(mswdr(end),'%.2f')])
    xlim([min(nr) max(nr)]);
    ylim([0 3])
    legend({'Weighted Mean','MSWD Test','Youngest Zircon','Bayesian estimate','Bayesian (uniform prior)'},'Location','northwest')
    % formatfigure;

end

%% Plot mean relative error for each dt_sigma
figure('OuterPosition',[600,0,400,1024]);

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
    
    subplot(length(dts),1,i); hold on; plot(repmat(nr,1,5), datar(:,6:10));
    plot(nr,ones(size(nr)),'--k');
    set(gca,'Xscale','log')
    xlabel('Number of zircons')
    ylabel('Average error / estimated error')
    title(['Delta t  =  ' num2str(dts(i)) ' sigma, MSWD: ' num2str(mswdr(end),'%.2f')])
    xlim([min(nr) max(nr)]);
    ylim([0 6])
    % formatfigure;
end

%% Plot MSWD as a function of N for different dt_sigma
cd('~/Desktop/Bayesian cryst timescales/results/');

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

%% Plot example datasets at different dt/sigma
cd ~/Desktop/'Bayesian cryst timescales'/

% Melts zircon distribution
dist = [0.000250810042630289,0.209249326535399,0.418247843028168,0.627246359520935,0.836373869090645,1.04374357155004,1.24322922919823,1.43183209158837,1.61129844550338,1.77688486888092,1.91661079082960,2.02685847965723,2.09922776245012,2.14214951000388,2.16049597963371,2.15970102681196,2.14446347280944,2.11932089388414,2.08857186945567,2.04960073596920,2.01041078899474,1.96673182187778,1.92474422194058,1.87915814526374,1.83509405646190,1.79113931931073,1.74610988807671,1.70247028643145,1.66040197025845,1.61746239633646,1.57521575240718,1.53439000675662,1.49515411445647,1.45616189914587,1.41776712500444,1.37956742674343,1.34262162411021,1.30694333629830,1.27238055574041,1.23881668637306,1.20536150622499,1.17185978304267,1.13930634959523,1.10770181433460,1.07702403102922,1.04755288935075,1.01891846438717,0.991307020548356,0.964597559566676,0.937935447447592,0.912547254896304,0.888432981912822,0.864397309803580,0.841090957191169,0.818607007784804,0.796806321622935,0.775811959280348,0.755412283398949,0.735655522179200,0.716614802353539,0.698283859993860,0.680538233670485,0.663358513619274,0.646723719544709,0.630552743575445,0.615079934324004,0.599677772774101,0.585285935569870,0.570974188360921,0.557289930703817,0.544044974206469,0.531081757193001,0.518528309854683,0.506357798936355,0.494356584164449,0.482751168373747,0.471302277398808,0.460045602813051,0.449055447325410,0.438327785548115,0.427997270253091,0.417817958928498,0.407979196700429,0.398511598892507,0.389449671179037,0.380749356382449,0.372391900054805,0.364322271587003,0.356549119619268,0.349118200705518,0.341882681922588,0.334822491050461,0.327796643262725,0.320830884764955,0.313845271864816,0.307017155215140,0.300244336445030,0.293473757055415,0.286703177665800,0.279932598276178];
xi = [0,0.0101010101010101,0.0202020202020202,0.0303030303030303,0.0404040404040404,0.0505050505050505,0.0606060606060606,0.0707070707070707,0.0808080808080808,0.0909090909090909,0.101010101010101,0.111111111111111,0.121212121212121,0.131313131313131,0.141414141414141,0.151515151515152,0.161616161616162,0.171717171717172,0.181818181818182,0.191919191919192,0.202020202020202,0.212121212121212,0.222222222222222,0.232323232323232,0.242424242424242,0.252525252525253,0.262626262626263,0.272727272727273,0.282828282828283,0.292929292929293,0.303030303030303,0.313131313131313,0.323232323232323,0.333333333333333,0.343434343434343,0.353535353535354,0.363636363636364,0.373737373737374,0.383838383838384,0.393939393939394,0.404040404040404,0.414141414141414,0.424242424242424,0.434343434343434,0.444444444444444,0.454545454545455,0.464646464646465,0.474747474747475,0.484848484848485,0.494949494949495,0.505050505050505,0.515151515151515,0.525252525252525,0.535353535353535,0.545454545454545,0.555555555555556,0.565656565656566,0.575757575757576,0.585858585858586,0.595959595959596,0.606060606060606,0.616161616161616,0.626262626262626,0.636363636363636,0.646464646464647,0.656565656565657,0.666666666666667,0.676767676767677,0.686868686868687,0.696969696969697,0.707070707070707,0.717171717171717,0.727272727272727,0.737373737373737,0.747474747474748,0.757575757575758,0.767676767676768,0.777777777777778,0.787878787878788,0.797979797979798,0.808080808080808,0.818181818181818,0.828282828282828,0.838383838383838,0.848484848484849,0.858585858585859,0.868686868686869,0.878787878787879,0.888888888888889,0.898989898989899,0.909090909090909,0.919191919191919,0.929292929292929,0.939393939393939,0.949494949494950,0.959595959595960,0.969696969696970,0.979797979797980,0.989898989898990,1];

dts = [0.1 0.5 1 2 5 10];

nzircs = 10;
agemax = 100.5;
agemin = 100;
agerange = agemax - agemin;

for i=1:length(dts)
    dt_sigma =dts(i);
    
    ages = agerange.*(1-randsample(xi,nzircs,true,dist)) + agemin;
    ages = sort(ages)';

    uncert = agerange./dt_sigma.*ones(size(ages));
    observedages = ages + randn(size(ages)).*uncert;
    
    sorted = sortrows([observedages, uncert],1);
    observedages = sorted(:,1);
    observeduncert = sorted(:,2);
    
    figure; hold on;
    errorbar(observedages,2*observeduncert,'.')
    plot([0,nzircs+1],[agemin agemin]);
    plot([0,nzircs+1],[agemax agemax]);

    [~,~,mswd]=wmean(observedages,observeduncert);
    xlabel('');
    xlim([0,nzircs+1]);
    set(gca,'Xtick',[]);
    title(['Delta t  =  ' num2str(dts(i)) ' sigma, MSWD = ' num2str(mswd,2)])
    % formatfigure;
end


%% Plot example datasets at different N
cd ~/Desktop/'Bayesian cryst timescales'/
if ~exist('dist','var'); load zircondistribution; end
Ns = [1 10 100];
dt_sigma = 1;

agemax = 100.1;
agemin = 100;
agerange = agemax - agemin;

for i=1:length(Ns)
    nzircs = Ns(i);
    
    ages = agerange.*(1-randsample(xi,nzircs,true,dist)) + agemin;
    ages = sort(ages)';

    uncert = agerange./dt_sigma.*ones(size(ages));
    observedages = ages + randn(size(ages)).*uncert;
    
    sorted = sortrows([observedages, uncert],1);
    observedages = sorted(:,1);
    observeduncert = sorted(:,2);
    
    figure; hold on;
    errorbar(observedages,2*observeduncert,'.')
    plot([0,nzircs+1],[agemin agemin]);
    plot([0,nzircs+1],[agemax agemax]);
    xlabel('');
    xlim([0,nzircs+1]);
    set(gca,'Xtick',[]);
    title(['Delta t  =  ' num2str(dt_sigma) ' sigma, N = ' num2str(nzircs)])
    % formatfigure;
end