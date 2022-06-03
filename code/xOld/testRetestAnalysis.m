%% testRetestAnalysis
% This script loads a blink data set into a MATLAB table variable. When
% run, it will aggregate data for a given subject and parameter(s), split
% by session. It will then produce a plot which describes the within session
% correlation for each session and will calculate the test retest reliability
% between sessions.
%% set up parameters

% load file path
dataPath = fileparts(fileparts(mfilename('fullpath')));
spreadsheet ='UPenn Ipsi Summary_25ms_02062022.csv';

% select blink features to analyze
varNamesToPlot = {'aucI', 'latencyI', 'timeUnderI', 'openTimeI', 'initVelocityI', ...
     'closeTimeI', 'maxClosingVelocityI', 'maxOpeningVelocityI', 'blinkRate'};

% decide whether analysis should include all 5 PSI levels or just the top
% three and select the appropriate subject list
highestOnly = true;
if highestOnly
    subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    14590, 14589, 14588, 14587, 14586};
else
    subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    14587, 14586};
end

xFit = linspace(log10(3),log10(70),50);

% create MATLAB table variable
T = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames = T.Properties.VariableNames;

%% session 1 and 2 subject blink feature across puff pressure
for vv = 1:length(varNamesToPlot)
    
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
        scans = scans(ismember(scans.numIpsi,(3:8)),:);

        % separate scans into a table for each of the sessions
        dates = unique(scans.scanDate);
        if highestOnly
           A = scans(ismember(scans.intendedPSI, 15),:);
           B = scans(ismember(scans.intendedPSI, 30),:);
           C = scans(ismember(scans.intendedPSI, 60),:);
           scans = vertcat(A, B, C);
        end
        sessOne = scans(ismember(scans.scanDate,dates(1,1)),:);
        sessTwo = scans(ismember(scans.scanDate,dates(2,1)),:);
        ii = find(strcmp(varNamesToPlot{vv},allVarNames));

        % session one data
        plotNum = plotNum + 1;
        y = sessOne.(allVarNames{ii});
        goodPoints = ~isnan(y);
        x = log10(sessOne.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        weights = sessOne.numIpsi;
        mSize = weights*20;
        
        % session one means
        sm = NaN(1,5);
        aa = sessOne(ismember(sessOne.intendedPSI, 3.5),:);
        bb = sessOne(ismember(sessOne.intendedPSI, 7.5),:);
        cc = sessOne(ismember(sessOne.intendedPSI, 15),:);
        dd = sessOne(ismember(sessOne.intendedPSI, 30),:);
        ee = sessOne(ismember(sessOne.intendedPSI, 60),:);
        sm(1) = mean(aa.(allVarNames{ii}), 'omitnan');
        sm(2) = mean(bb.(allVarNames{ii}), 'omitnan');
        sm(3) = mean(cc.(allVarNames{ii}), 'omitnan');
        sm(4) = mean(dd.(allVarNames{ii}), 'omitnan');
        sm(5) = mean(ee.(allVarNames{ii}), 'omitnan');

        % make plot
        subplot(2,length(subList),plotNum);
        scatter(x,y,mSize);
%         hold on
%         scatter(log10([3.5 7.5 15 30 60]),sm,300);
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
        hold on
        plot(x,fitObj.Fitted,'-r');
        xlim(log10([2 100]));
        pX(end+1) = fitObj.Coefficients.Estimate(2);
        oX(end+1) = fitObj.Coefficients.Estimate(1);
        rsquare = fitObj.Rsquared.Ordinary;
        if rsquare > 1 || rsquare < 0
            rsquare = nan;
        end
        title(['Subject ' num2str(subList{ss})], 'FontSize', 14)
        if plotNum ~= 1
            yticklabels("");
            xticklabels("");
            xticks([]);
            yticks([]);
        else
            ylabel(['Session one ' varNamesToPlot{vv}], 'FontSize', 14)
            xlabel('puff pressure [log psi]', 'FontSize', 14)
        end

        % session two data
        y = sessTwo.(allVarNames{ii});
        goodPoints = ~isnan(y);
        x = log10(sessTwo.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        weights = sessTwo.numIpsi;
        mSize = weights*20;

        % make plot
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
        if plotNum == 1
            ylabel(['Session two ' varNamesToPlot{vv}], 'FontSize', 14)
        end
        yticklabels("");
        xticklabels("");
        xticks([]);
        yticks([]);
    end
    
    % plot parameter test retest values across subjects
    figure();
    pl = subplot(1,1,1);
    plot(pX, pY, 'ob', 'MarkerSize', 10);
    fitObj = fitlm(pX,pY,'RobustOpts', 'on');
    hold on
%     plot(pX,fitObj.Fitted,'-r')
    plot(fitObj);
%     plot((-2:5),(-2:5),'k');
    pl.Box = 'off';
    rsquare = fitObj.Rsquared.Ordinary;
    if rsquare > 1 || rsquare < 0
        rsquare = nan;
    end
    title([varNamesToPlot{vv} ' slope by session - ' sprintf(' R^2=%2.2f',rsquare)], 'FontSize', 16)
    xlabel(['Slope'], 'FontSize', 16)
    ylabel(['Slope'], 'FontSize', 16)
    ylim(xlim);
    axis(pl, 'square');
    
    % plot offset test retest values across subjects
    figure();
    pl = subplot(1,1,1);
    plot(oX, oY, 'ob', 'MarkerSize', 10);
    fitObj = fitlm(oX,oY,'RobustOpts', 'on');
%     hold on
%     plot(oX,fitObj.Fitted,'-r')
    plot(fitObj);
%     plot((-2:14),(-2:14),'k');
    pl.Box = 'off';
    rsquare = fitObj.Rsquared.Ordinary;
    if rsquare > 1 || rsquare < 0
        rsquare = nan;
    end
    title([varNamesToPlot{vv} ' offset by session - ' sprintf(' R^2=%2.2f',rsquare)], 'FontSize', 16)
    xlabel(['Offset'], 'FontSize', 16)
    ylabel(['Offset'], 'FontSize', 16)
    ylim(xlim);
    axis(pl, 'square');
    
end

%% test retest slope or offset Bland Altman analysis
figure();

for vv = 1:length(varNamesToPlot)
    
    offset = false;
    
    oX = [];
    oY = [];
    pX = [];
    pY = [];
    
    for ss = 1:length(subList)

        % find scans for desired subject
        scans = T(ismember(T.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);
        scans = scans(ismember(scans.numIpsi,(3:8)),:);

        % separate scans into a table for each of the sessions
        dates = unique(scans.scanDate);
        if highestOnly
           A = scans(ismember(scans.intendedPSI, 15),:);
           B = scans(ismember(scans.intendedPSI, 30),:);
           C = scans(ismember(scans.intendedPSI, 60),:);
           scans = vertcat(A, B, C);
        end
        sessOne = scans(ismember(scans.scanDate,dates(1,1)),:);
        sessTwo = scans(ismember(scans.scanDate,dates(2,1)),:);
        ii = find(strcmp(varNamesToPlot{vv},allVarNames));

        % session one data
        y = sessOne.(allVarNames{ii});
        goodPoints = ~isnan(y);
        x = log10(sessOne.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        weights = sessOne.numIpsi;
        mSize = weights*20;
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
        oX(end+1) = fitObj.Coefficients.Estimate(1);
        pX(end+1) = fitObj.Coefficients.Estimate(2);

        % session two data
        y = sessTwo.(allVarNames{ii});
        goodPoints = ~isnan(y);
        x = log10(sessTwo.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        weights = sessTwo.numIpsi;
        mSize = weights*20;
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
        oY(end+1) = fitObj.Coefficients.Estimate(1);
        pY(end+1) = fitObj.Coefficients.Estimate(2);
    end
    
    % Bland-Altman Analysis
    if offset
        
        % get BA stats
        differences = oX - oY;
        meanD = mean(differences,'omitnan');
        meanO = (oX + oY) ./ 2;
        
        % get IQR
        range = iqr(differences);
        low = meanD - 0.5*range;
        high = meanD + 0.5*range;
        
        % plot stats
        scatter(meanO, differences, 150, 'MarkerEdgeColor', 'r');
        hold on
        yline(meanD,'-',num2str(meanD));
        yline(low,'--',num2str(low));
        yline(high,'--',num2str(high));
        title([varNamesToPlot{vv}], 'FontSize', 16)
        xlabel(['Mean offset'], 'FontSize', 16)
        ylabel(['Offset difference'], 'FontSize', 16)
        axis(pl, 'square');
    else
        
        % get BA stats
        differences = pX - pY;
        meanD = mean(differences,'omitnan');
        meanO = (pX + pY) ./ 2;
        
        % get IQR
        range = iqr(differences);
        low = meanD - 0.5*range;
        high = meanD + 0.5*range;
        
        % plot stats
        scatter(meanO, differences, 150, 'MarkerEdgeColor', 'r');
        hold on
        yline(meanD,'-',num2str(meanD));
        yline(low,'--',num2str(low));
        yline(high,'--',num2str(high));
        title([varNamesToPlot{vv}], 'FontSize', 16)
        xlabel(['Mean slope'], 'FontSize', 16)
        ylabel(['Slope difference'], 'FontSize', 16)
        axis(pl, 'square');
    end
    
end

%% slope vs offset

figure();
feature = 'maxClosingVelocityI';
    
pY = [];
oX = [];

for ss = 1:length(subList)

    % find scans for desired subject
    scans = T(ismember(T.subjectID,subList{ss}),:);
    scans = scans(ismember(scans.valid,'TRUE'),:);
    scans = scans(ismember(scans.numIpsi,(3:8)),:);

    % separate scans into a table for each of the sessions
    if highestOnly
       A = scans(ismember(scans.intendedPSI, 15),:);
       B = scans(ismember(scans.intendedPSI, 30),:);
       C = scans(ismember(scans.intendedPSI, 60),:);
       scans = vertcat(A, B, C);
    end
    ii = find(strcmp(varNamesToPlot{vv},allVarNames));

    % subject data
    y = scans.(allVarNames{ii});
    goodPoints = ~isnan(y);
    x = log10(scans.PSI);
    x = x(goodPoints);
    y = y(goodPoints);
    [x,idxX]=sort(x);
    y = y(idxX);
    weights = scans.numIpsi;
    mSize = weights*20;
    fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
    oX(end+1) = fitObj.Coefficients.Estimate(2)*median(x)+fitObj.Coefficients.Estimate(1);
    pY(end+1) = fitObj.Coefficients.Estimate(2);

end

% plot parameter test retest values across subjects
pl = subplot(1,1,1);
plot(oX, pY, 'ob', 'MarkerSize', 10);
fitObj = fitlm(oX,pY,'RobustOpts', 'on');
hold on
plot(fitObj);
% plot((-2:5),(-2:5),'k');
pl.Box = 'off';
rsquare = fitObj.Rsquared.Ordinary;
if rsquare > 1 || rsquare < 0
    rsquare = nan;
end
title([feature ' slope vs. offset ' sprintf(' R^2=%2.2f',rsquare)], 'FontSize', 16)
xlabel(['Offset'], 'FontSize', 16)
ylabel(['Slope'], 'FontSize', 16)
axis(pl, 'square');