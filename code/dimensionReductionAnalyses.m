%% dimensionReductionAnalyses
% This script loads a blink data set into a MATLAB table variable. When
% run, it will aggregate data for a given subject for the given parameters
% across sessions. It will create a matrix of slopes for each parameter.
clear all

% load file path
dataPath = fileparts(fileparts(mfilename('fullpath')));
spreadsheet ='UPenn Ipsi Summary_25ms_02062022.csv';

% choose subject and parameters
highestOnly = true;
if highestOnly
    subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
        14590, 14589, 14588, 14587, 14586};
else
    subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
        14587, 14586};
end
varNamesToPlot = {'auc', 'latency', 'timeUnder', 'openTime', 'initVelocity', ...
    'closeTime', 'maxClosingVelocity', 'maxOpeningVelocity', 'blinkRate'};

varIdxToUse = 1:9;

% create MATLAB table variable
T = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames = T.Properties.VariableNames;
slopes = zeros(length(subList), length(varNamesToPlot));
slopesSessOne = zeros(length(subList), length(varNamesToPlot));
slopesSessTwo = zeros(length(subList), length(varNamesToPlot));
offsets = zeros(length(subList), length(varNamesToPlot));
offsetsSessOne = zeros(length(subList), length(varNamesToPlot));
offsetsSessTwo = zeros(length(subList), length(varNamesToPlot));

%% create slopes matrix containing the slope values for each var and subject
for vv = 1:length(varNamesToPlot)

    for ss = 1:length(subList)

        % find scans for desired subject
        scans = T(ismember(T.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);
        scans = scans(ismember(scans.numIpsi,(3:8)),:);
        if highestOnly
            A = scans(ismember(scans.intendedPSI, 15),:);
            B = scans(ismember(scans.intendedPSI, 30),:);
            C = scans(ismember(scans.intendedPSI, 60),:);
            scans = vertcat(A, B, C);
        end
        ii = find(strcmp(varNamesToPlot{vv},allVarNames));
        weights = scans.numIpsi;
        dates = unique(scans.scanDate);
        sessOne = scans(ismember(scans.scanDate,dates(1,1)),:);
        sessTwo = scans(ismember(scans.scanDate,dates(2,1)),:);
        weightsSessOne = sessOne.numIpsi;
        weightsSessTwo = sessTwo.numIpsi;

        % subject parameter data across sessions
        y = scans.(allVarNames{ii});
        goodPoints = ~isnan(y);
        x = log10(scans.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
        slopes(ss, vv) = fitObj.Coefficients.Estimate(2);

        % Get the y value at median x and it will be our offset
        offsets(ss,vv) = fitObj.Coefficients.Estimate(2)*median(x)+fitObj.Coefficients.Estimate(1);

        % subject parameter data session 1
        y = sessOne.(allVarNames{ii});
        goodPoints = ~isnan(y);
        x = log10(scans.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weightsSessOne);
        slopesSessOne(ss, vv) = fitObj.Coefficients.Estimate(2);

        % Get the y value at median x and it will be our offset
        offsetsSessOne(ss,vv) = fitObj.Coefficients.Estimate(2)*median(x)+fitObj.Coefficients.Estimate(1);

        % subject parameter data session 2
        y = sessTwo.(allVarNames{ii});
        goodPoints = ~isnan(y);
        x = log10(scans.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weightsSessTwo);
        slopesSessTwo(ss, vv) = fitObj.Coefficients.Estimate(2);

        % Get the y value at median x and it will be our offset
        offsetsSessTwo(ss, vv) = fitObj.Coefficients.Estimate(2)*median(x)+fitObj.Coefficients.Estimate(1);
    end

end

comp = corr(slopes);

%% Make test-retest analysis.

% First we do 2 PCA analyses on the slope and offset session aggregate
% variables. We then project session 1 and session 2 vectors onto the PC1
% and PC2 vectors and correlate their new scores.

% Save measurements and names in cells
allMeasures = {slopes, offsets; ...
    slopesSessOne, slopesSessTwo; ...
    offsetsSessOne, offsetsSessTwo};
allMeasureNames = {'Slopes', 'Offsets', ...
    'slopesSessOne', 'slopesSessTwo', ...
    'offsetsSessOne', 'offsetsSessTwo'};

% Initialize figures
figure1 = figure('Renderer', 'painters', 'Position', [164 71 1401 891]);
figure2 = figure('Renderer', 'painters', 'Position', [164 71 1401 891]);
plotCounter = 1;
plotCounter2 = 1;
combinedSessions = [];

% Loop over the set [slopes, offsets, (slopes offsets)], and for each data set perform an NMF
% decomposition. Save the nmfResults struct and plot some diagnostics

for ii = 1:2

    % Standardize the measures
    standardized = (allMeasures{1,ii}-mean(allMeasures{1,ii}))./std(allMeasures{1,ii});
    varNames = strcat(varNamesToPlot(varIdxToUse),allMeasureNames{ii});

    % Conduct an NMF dimension reduction for the aggregate measurements.
    % First, loop over the number of dimensions and save the RMS residual
    nmfResidual(1) = nan;
    for dd = 2:4
        [~,~,nmfResidual(dd)] = nnmf(standardized(:,varIdxToUse),dd);
    end

    % We end up running this with 3 dimensions
    [W,H] = nnmf(standardized(:,varIdxToUse),3);
    coeff = H';
    score = W;

    % Save NMF results to a structure
    nmfResults.(allMeasureNames{1,ii}).('coeff') = coeff;
    nmfResults.(allMeasureNames{1,ii}).('score') = score;
    nmfResults.(allMeasureNames{1,ii}).('standardized') = standardized;

    % Create some diagnostic plots

    % Residual by dimensions
    set(0,'CurrentFigure',figure1)
    subplot(2, 3, plotCounter)
    plot(1:4,nmfResidual)
    xlabel('number of dimensions'); ylabel('RMS residual');

    % biplot
    subplot(2, 3, plotCounter+1)
    biplot(coeff(:,1:3),'scores',score(:,1:3),'varLabels',varNames)
    title(allMeasureNames{1,ii})
    axis equal

    % scores
    subplot(2, 3, plotCounter+2)
    scatter(score(:,1),score(:,2),'or')
    axis square
    title(allMeasureNames{1,ii})
    xlabel('component 1'); ylabel('component 2'); zlabel('component 3');

    % Project the individual sessions to the NMF solution
    sessionOneStandard = (allMeasures{ii+1,1}-mean(allMeasures{ii+1,1}))./std(allMeasures{ii+1,1});
    sessionTwoStandard = (allMeasures{ii+1,2}-mean(allMeasures{ii+1,2}))./std(allMeasures{ii+1,2});
    sessionOneProjected = sessionOneStandard(:,varIdxToUse)*coeff(:,1:2);
    sessionTwoProjected = sessionTwoStandard(:,varIdxToUse)*coeff(:,1:2);

    set(0,'CurrentFigure',figure2)
    subplot(2,2,plotCounter2)
    mdl = fitlm(sessionOneProjected(:,1), sessionTwoProjected(:,1));
    plot(mdl, 'Marker', 'o', 'MarkerEdgeColor','b', 'MarkerFaceColor','b')
    legend off
    xlabel('Session 1 projected onto NMF1')
    ylabel('Session 2 projected onto NMF1')
    title(allMeasureNames{ii})
    axis equal
    xlim([-4 4])
    ylim([-4 4])
    axis square

    subplot(2,2,plotCounter2+1)
    mdl = fitlm(sessionOneProjected(:,2), sessionTwoProjected(:,2));
    plot(mdl, 'Marker', 'o', 'MarkerEdgeColor','b', 'MarkerFaceColor','b')
    legend off
    xlabel('Session 1 projected onto NMF2')
    ylabel('Session 2 projected onto NMF2')
    title(allMeasureNames{ii})
    axis equal
    xlim([-4 4])
    ylim([-4 4])
    axis square

    plotCounter = plotCounter + 3;
    plotCounter2 = plotCounter2 + 2;
end

%% Correlation of NMF dimensions derived from slopes and offsets
slopesNMF1 = nmfResults.Slopes.standardized(:,varIdxToUse) * nmfResults.Slopes.coeff(:,1);
slopesNMF2 = nmfResults.Slopes.standardized(:,varIdxToUse) * nmfResults.Slopes.coeff(:,2);
offsetNMF1 = nmfResults.Offsets.standardized(:,varIdxToUse) * nmfResults.Offsets.coeff(:,1);
offsetNMF2 = nmfResults.Offsets.standardized(:,varIdxToUse) * nmfResults.Offsets.coeff(:,2);
