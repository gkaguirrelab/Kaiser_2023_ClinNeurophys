%% compareParams
% This script loads a blink data set into a MATLAB table variable. When
% run, it will aggregate data for a given subject for 2 parameters across
% sessions. It will then produce a plot which describes the across session
% trend for each subject and each of the parameters. It will then produce
% a second plot which describes the correlation of the slope for the two 
% parameters across subjects.
%%

% load file path
dataPath = fileparts(fileparts(mfilename('fullpath')));
spreadsheet ='2_2022.csv';

% choose subject and parameters
subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591};
varNamesToPlot = {'aucI', 'latencyI'};

xFit = linspace(log10(3),log10(70),50);
ylims = {[0 5e4], [30 80]};

% create MATLAB table variable
T = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames = T.Properties.VariableNames;

figure();
plotNum = 0;
pX = [];
pY = [];
oX = [];
oY = [];

for ss = 1:length(subList)

    % find scans for desired subject
    scans = T(ismember(T.subjectID,subList{ss}),:);
    scans = scans(ismember(scans.valid,'TRUE'),:);
    ii = find(strcmp(varNamesToPlot{1},allVarNames));
    pp = find(strcmp(varNamesToPlot{2},allVarNames));
    weights = scans.numIpsi;
    mSize = weights*20;

    % subject parameter one data
    plotNum = plotNum + 1;
    y = scans.(allVarNames{ii});
    goodPoints = ~isnan(y);
    x = log10(scans.PSI);
    x = x(goodPoints);
    y = y(goodPoints);
    [x,idxX]=sort(x);
    y = y(idxX);

    % make paramater one plot
    subplot(2,length(subList),plotNum);
    scatter(x,y,mSize);
    fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
    hold on
    plot(x,fitObj.Fitted,'-r')
    xlim(log10([2 100]));
    pX(end+1) = fitObj.Coefficients.Estimate(2);
    oX(end+1) = fitObj.Coefficients.Estimate(1);
    rsquare = fitObj.Rsquared.Ordinary;
    if rsquare > 1 || rsquare < 0
        rsquare = nan;
    end
    title([varNamesToPlot{1} ' - ' num2str(subList{ss}) sprintf(' R^2=%2.2f',rsquare)])
    xlabel('puff pressure [log psi]')
    ylim(ylims{1});

    % subject parameter two data
    y = scans.(allVarNames{pp});
    goodPoints = ~isnan(y);
    x = log10(scans.PSI);
    x = x(goodPoints);
    y = y(goodPoints);
    [x,idxX]=sort(x);
    y = y(idxX);
    
    % make paramater two plot
    subplot(2,length(subList),plotNum + length(subList));
    scatter(x,y,mSize);
    fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
    hold on
    plot(x,fitObj.Fitted,'-r')
    xlim(log10([2 100]));
    pY(end+1) = fitObj.Coefficients.Estimate(2);
    oY(end+1) = fitObj.Coefficients.Estimate(1);
    rsquare = fitObj.Rsquared.Ordinary;
    if rsquare > 1 || rsquare < 0
        rsquare = nan;
    end
    title([varNamesToPlot{2} ' - ' num2str(subList{ss}) sprintf(' R^2=%2.2f',rsquare)])
    xlabel('puff pressure [log psi]')
    ylim(ylims{2});

end

% plot parameter test retest values across subjects
figure();
subplot(1,2,1);
plot(pX,pY,'ob');
fitObj = fitlm(pX,pY,'RobustOpts', 'on');
hold on
plot(pX,fitObj.Fitted,'-r')
rsquare = fitObj.Rsquared.Ordinary;
if rsquare > 1 || rsquare < 0
    rsquare = nan;
end
title(['Across subject parameter - ' sprintf(' R^2=%2.2f',rsquare)])
xlabel([varNamesToPlot{1} ' parameter'])
ylabel([varNamesToPlot{2} ' parameter'])

% plot offset test retest values across subjects
subplot(1,2,2);
plot(oX,oY,'ob');
fitObj = fitlm(oX,oY,'RobustOpts', 'on');
hold on
plot(oX,fitObj.Fitted,'-r')
rsquare = fitObj.Rsquared.Ordinary;
if rsquare > 1 || rsquare < 0
    rsquare = nan;
end
title(['Across subject offset - ' sprintf(' R^2=%2.2f',rsquare)])
xlabel([varNamesToPlot{1} ' offset'])
ylabel([varNamesToPlot{2} ' offset'])