% Ensure you are running from the right folder!
% The code must be able to find the CSV files for each sample..

%% Run Metropolis walker to find eruption age of each sample
nsteps = 500000;
burnin = 10000;

% Allocate struct arrays for output
sample.elements = {'Name','Distribution','LiteratureAge','EruptionAge','EruptionAge_Sigma','SaturationAge','SaturationAge_Sigma'}';


sample.Name = {'CrowleyBishopTuff';'CrowleyFilteredBishopTuff';'IckertBishopTuff';'WotzlawFCTMLX';'WotzlawFCTMLXYoung';'WotzlawNutrasCreekDacite';'WotzlawFishCanyonAll';'WotzlawKilgoreTuff';'RiveraAlderCreek';};
for i=2:length(sample.elements)
    sample.(sample.elements{i}) = NaN(size(sample.Name));
end
sample.ReportedWMAge = [0.7679; 0.7671; 0.7671; NaN; NaN; NaN; NaN; 4.4876; 1.1978;]; 
sample.ReportedWMUncertainty = [0.000725; 0.0009/2; 0.0009/2; NaN; NaN; NaN; NaN; 0.0023; 0.0046/2;]; 
sample.ReportedArAge = [0.7666; 0.7666; 0.7666; 28.201; 28.201; 28.201; 28.201; NaN; 1.1848;]; 
sample.ReportedArUncertainty = [0.0008/2; 0.0008/2; 0.0008/2; 0.046/2; 0.046/2; 0.046/2; 0.046/2; NaN; 0.0006;]; 
sample.Distribution = {'VolcanicZirconLowX';'VolcanicZirconLowX';'VolcanicZirconLowX';'VolcanicZircon';'VolcanicZircon';'VolcanicZircon';'VolcanicZircon';'VolcanicZircon';'VolcanicZircon';};
% sample.Distribution = {'Uniform';'Uniform';'Uniform';'Uniform';'Uniform';'Uniform';'Uniform';'Uniform';'Uniform';};

% Run Metropolis walker for each sample
for i=1:length(sample.Name)
    %%%%%%%%%%%%%%%%%%%%%%%%% Bootstrapped prior %%%%%%%%%%%%%%%%%%%%%%%%%%

    % Load raw ages
        data = load([sample.Name{i} '.csv']);
        sample.(sample.Name{i}).Age = data(:,1);
        sample.(sample.Name{i}).Age_Sigma = data(:,2);
        agescaled = sample.(sample.Name{i}).Age - min(sample.(sample.Name{i}).Age);
        agescaled = agescaled./max(agescaled);
    % Obtain ksdensity values between 0 and 1+bandwidth
        npoints = 100;
        bwcutoff = 2;
        [f,x,bw] = ksdensity(agescaled);
        BootstrappedDistribution = interp1(x,f,linspace(0,1+bwcutoff*bw,npoints));
    % Normalize the distribution
        x = linspace(0,1,npoints);
        BootstrappedDistribution = BootstrappedDistribution./trapz(x,BootstrappedDistribution);
        exportmatrix(flipud(BootstrappedDistribution'),'BootstrappedDistribution.tsv','\t');
        % figure; plot(x,BootstrappedDistribution)
        % set(gca,'xdir','reverse')
    system(['../../src/tzircrystmetropolis ' num2str(nsteps) ' ./BootstrappedDistribution.tsv ' sample.Name{i} '.csv > metropolisdata.tsv']);

    %%%%%%%%%%%%%%%%%%%%%% Uniform or Volcanic prior %%%%%%%%%%%%%%%%%%%%%%

    system(['../../src/tzircrystmetropolis ' num2str(nsteps) ' ../../distributions/' sample.Distribution{i} 'Distribution.tsv ' sample.Name{i} '.csv > metropolisdata.tsv']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    load metropolisdata.tsv
    % figure; plot(metropolisdata(:,1:2))

    sample.SaturationAge(i) = nanmean(metropolisdata(burnin:end,2));
    sample.SaturationAge_Sigma(i) = nanstd(metropolisdata(burnin:end,2));
    sample.EruptionAge(i) = nanmean(metropolisdata(burnin:end,1));
    sample.EruptionAge_Sigma(i) = nanstd(metropolisdata(burnin:end,1));
end
exportdataset(sample,'BayesianResults.csv',',');


%% Rank-order plots for each

c = lines(5);
for i=1:length(sample.Name)
    data = load([sample.Name{i} '.csv']);
    data = sortrows(data);
    figure; errorbar(1:size(data,1),data(:,1),data(:,2),'.','Color',c(1,:))
    hold on; plot([0, size(data,1)+1],[sample.SaturationAge(i),sample.SaturationAge(i)],'Color',c(2,:))
    hold on; plot([0, size(data,1)+1],[sample.EruptionAge(i),sample.EruptionAge(i)],'Color',c(3,:))
    hold on; errorbar(size(data,1)+1,sample.EruptionAge(i),sample.EruptionAge_Sigma(i)*2,'.','Color',c(3,:))
    hold on; errorbar(size(data,1)+2,sample.ReportedWMAge(i),sample.ReportedWMUncertainty(i)*2,'.','Color',c(4,:))
    hold on; errorbar(size(data,1)+3,sample.ReportedArAge(i),sample.ReportedArUncertainty(i)*2,'.','Color',c(5,:))
    ylabel('Age (Ma)')
    ymin = min([data(:,1); sample.EruptionAge(i)]);
    ymax = max([data(:,1); sample.SaturationAge(i)]);
    xlim([0 size(data,1)+4])
    title(regexprep(sample.Name{i},'_','-'));
    formatfigure;
%     saveas(gcf,[sample.Name{i} 'ages.eps'],'epsc');
end

