function LaminarPlot()

f = figure('color','w');

%===============================================================================
ax = subplot(1,3,1);
[layersExamined,corrs,pVals,pValsCorr,cellTypes] = TheBestLayers('cell',false);
[isSig_i,isSig_j] = find(pValsCorr < 0.05);
numSig = length(isSig_i);
for k = 1:numSig
    ii = isSig_i(k);
    jj = isSig_j(k);
    fprintf(1,'%s (%s):r = %.2f, p_corr = %.3g\n',cellTypes{jj},layersExamined{ii},corrs(ii,jj),pValsCorr(ii,jj));
end
LabelCurrentAxes('A',ax,18,'topLeft');
colorbar('off')

%===============================================================================
ax = subplot(1,3,2:3);
TheBestLayers('genes',false);
LabelCurrentAxes('B',ax,18,'topLeft');
f.Position = [626,875,1296,328];

end
