% Run from /examples/BootstrappedSaturationDistributionExamples
%% Load and stack Bergell data
% Data from Samperton et al., 2017
bergell = str2double(importc('BergellAll.csv',','));
% Plot data
figure; plot(bergell(:),'.','MarkerSize',10)
xlabel("N"); ylabel("Age, Ma")

for i=1:size(bergell,2)
    bergell(:,i) = bergell(:,i) - min(bergell(:,i));
    bergell(:,i) = bergell(:,i)./max(bergell(:,i));
end
% Plot scaled data
figure; plot(bergell(:),'.','MarkerSize',10)
xlabel("N"); ylabel("Age, scaled")

%% Load and stack Elba data
% load data -- from Barboni, Annen, and Schoene 2015, including many units
elba = str2double(importc('ElbaAll.csv',','));
% Plot data
figure; plot(elba(:),'.','MarkerSize',10)
xlabel("N"); ylabel("Age, Ma")
for i=1:size(elba,2)
    elba(:,i) = elba(:,i) - min(elba(:,i));
    elba(:,i) = elba(:,i)./max(elba(:,i));
end
% Plot scaled data
figure; plot(elba(:),'.','MarkerSize',10)
xlabel("N"); ylabel("Age, scaled")

%% Load and stack Fish canyon data
% Data from Wotzlaw et al. 2014 -- three units are Nutras Creek Dacit,  Fish Canyon Tuff, and andesite enclaves therein
FC = str2double(importc('FishCanyonAll.csv',','));
% Plot data
figure; plot(FC(:),'.','MarkerSize',10)

for i=1:size(FC,2)
    FC(:,i) = FC(:,i) - min(FC(:,i));
    FC(:,i) = FC(:,i)./max(FC(:,i));
end
% Plot scaled data
figure; plot(FC(:),'.','MarkerSize',10)
xlabel("N"); ylabel("Age, scaled")

%% Calculate and plot distributions
npoints = 100;
bwcutoff = 1;
[bergell_f,bergell_x,bergell_bw] = ksdensity(bergell_all);
[elba_f,elba_x,elba_bw] = ksdensity(elba_all);
[FC_f,FC_x,FC_bw] = ksdensity(FC_all);

% Obtain ksdensity values between 0-bandwith and 1+2*bandwidth
bdist = interp1(bergell_x,bergell_f,linspace(0-bwcutoff*bergell_bw,1+2*bergell_bw,npoints));
edist = interp1(elba_x,elba_f,linspace(0-bwcutoff*elba_bw,1+2*elba_bw,npoints));
FCdist = interp1(FC_x,FC_f,linspace(0-bwcutoff*FC_bw,1+2*FC_bw,npoints));

% Normalize each distribution
x = linspace(0,1,npoints);
bdist = bdist./trapz(x,bdist);
edist = edist./trapz(x,edist);
FCdist = FCdist./trapz(x,FCdist);

% Plot results, scaled by age
c = lines(5);
figure; plot(x,bdist,'Color',c(4,:))
hold on; plot(x,edist,'Color',c(5,:))
hold on; plot(x,FCdist,'Color',c(3,:))
legend('Samperton - Bergell','Barboni - Elba','Wotzlaw - Fish Canyon');

ylim([0 3]);
set(gca,'xdir','reverse')
xlabel("Age (scaled)")
ylabel("Density")
formatfigure

