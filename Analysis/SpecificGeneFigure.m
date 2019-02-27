function SpecificGeneFigure()

% Franklin-Paxinos atlas:
convertToFP = false;

% Regress out neuron density:
doNNRegression = false;

%-------------------------------------------------------------------------------
f = figure('color','w');
f.Position = [626,875,1296,328];

ax = subplot(1,3,1);
SpecificGenes(false,'XJ',false,true,convertToFP,doNNRegression);
LabelCurrentAxes('A',gca,18,'topLeft');
r = [-0.7,0.7];
xlim(r);caxis(r)
ax.Position = [0.0324    0.1100    0.3380    0.8150];

subplot(1,3,2)
SpecificGenes(false,'Pvalb',true,false,convertToFP,doNNRegression);
LabelCurrentAxes('B',gca,18,'topLeft');
legend('off'); title('')

subplot(1,3,3)
SpecificGenes(false,'Grin3a',true,false,convertToFP,doNNRegression);
LabelCurrentAxes('C',gca,18,'topLeft');
legend('off'); title('')

end
