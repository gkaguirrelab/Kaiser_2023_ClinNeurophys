dataPath = fileparts(fileparts(mfilename('fullpath')));

T = readtable(fullfile(dataPath,'data','gkaDemoBlinkCNSData.csv'));
varNames = T.Properties.VariableNames;

xFit = linspace(log10(4),log10(50),50);
    figure

for ii = 2:length(varNames)
    subplot(4,3,ii-1);
    x = log10(T.PSI);
    y = T.(varNames{ii});
    goodPoints = ~isnan(y);
    x = x(goodPoints);
    y = y(goodPoints);
    [x,idxX]=sort(x);
    y = y(idxX);
    plot(x,y,'ok');
    [fitObj,G] = L3P(x,y);
    hold on
    plot(xFit,fitObj(xFit),'-r')
    xlim(log10([4 50]));
    title(varNames{ii})
    xlabel('puff pressure [log psi]')
end