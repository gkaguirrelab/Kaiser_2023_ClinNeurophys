dataPath = fileparts(fileparts(mfilename('fullpath')));
fileName = 'gkaDemoBlinkCNSData.csv';
fileName = 'p100.csv';

T = readtable(fullfile(dataPath,'data',fileName));
varNames = T.Properties.VariableNames;

xFit = linspace(log10(4),log10(50),50);
    figure

for ii = 2:length(varNames)
    subplot(5,3,ii-1);
    x = log10(T.PSI);
    y = T.(varNames{ii});
    switch varNames{ii}
        case 'latency'
            y = 100-y;
            titleStr = '100-latency';
        otherwise
            titleStr = varNames{ii};
    end
    goodPoints = ~isnan(y);
    x = x(goodPoints);
    y = y(goodPoints);
    [x,idxX]=sort(x);
    y = y(idxX);
    plot(x,y,'ok');
    [fitObj,G] = L3P(x,y);
    hold on
    plot(xFit,fitObj(xFit),'-r')
    xlim(log10([2 70]));
    title([titleStr sprintf(' R^2=%2.2f',G.rsquare)])
    xlabel('puff pressure [log psi]')
end