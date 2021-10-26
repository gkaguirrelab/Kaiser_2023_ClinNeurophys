dataPath = fileparts(fileparts(mfilename('fullpath')));

subList = {'p96.csv','p97.csv','p99.csv','p100.csv'};

varNamesToPlot = {'latency','auc'};

xFit = linspace(log10(3),log10(70),50);

figure();

idx = 1;
for ss = 1:length(subList)
    T = readtable(fullfile(dataPath,'data',subList{ss}));
    allVarNames = T.Properties.VariableNames;
    
    for vv = 1:length(varNamesToPlot)
        ii = find(strcmp(varNamesToPlot{vv},allVarNames));
        
        y = T.(allVarNames{ii});
        goodPoints = ~isnan(y);
        x = log10(T.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        
        subplot(length(subList),length(varNamesToPlot),idx);
        idx = idx+1;
        
        plot(x,y,'ok');
        [fitObj,G] = L3P(x,y);
        hold on
        plot(xFit,fitObj(xFit),'-r')
        xlim(log10([2 100]));
        rsquare = G.rsquare;
        if rsquare > 1 || rsquare < 0
            rsquare = nan;
        end
        title([varNamesToPlot{vv} ' - ' subList{ss} sprintf(' R^2=%2.2f',rsquare)])
        xlabel('puff pressure [log psi]')
    end
end