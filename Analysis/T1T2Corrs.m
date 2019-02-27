function T1T2Corrs()

removeTitles = true;

f = figure('color','w');
f.Position = [626,875,1296,328];

subplot(141)
PlotCorrWithX('cytoarchitecture','C',removeTitles)
subplot(142)
PlotCorrWithX('PV_SST_L23','D',removeTitles)
subplot(143)
PlotCorrWithX('inStrength','E',removeTitles)
subplot(144)
PlotCorrWithX('HarrisHierarchy','F',removeTitles)

function PlotCorrWithX(whatFeature,whatLabel,removeTitles)
    MyelinCorrs(whatFeature,false)
    LabelCurrentAxes(whatLabel,gca,18,'topLeft');
    if removeTitles
        title('');
    end
end

end
