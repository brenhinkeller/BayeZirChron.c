function pt = getTransitionProbability(stepfactor)
workdir=['out' num2str(stepfactor)];
system(['mkdir -p ' workdir]);
system(['cd ' workdir '; ../tzircrystmetropolis 100 100000 ' num2str(stepfactor) ' ../MeltsTZircDistribtuion.csv ../zircondata.tsv > metropolisdata.tsv']);
load([workdir '/metropolisdata.tsv']);
pt = sum(diff(metropolisdata(:,3))~=0)./length(metropolisdata(:,3));