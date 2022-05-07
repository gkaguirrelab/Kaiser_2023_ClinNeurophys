%% dimensionReductionAnalyses
% This script loads a blink data set into a MATLAB table variable. When
% run, it will aggregate data for a given subject for the given parameters
% across sessions. It will create a matrix of slopes for each parameter.

% Housekeeping
clear; close all
rng('default')
warnState = warning();
warning('off','stats:nnmf:LowRank');
warning('off','stats:statrobustfit:IterationLimit');

% load file path
dataPath = fileparts(fileparts(mfilename('fullpath')));
spreadsheet ='UPenn Ipsi Summary_25ms_02062022.csv';

% choose subject and parameters
highestOnly = true;
if highestOnly
%     subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
%         14590, 14589, 14588, 14587, 14586};
    subList = {15512, 15507, 15506, 15505, 14595, 14594, 14593, 14592, 14591, ...
        14590, 14589, 14588, 14587, 14586};
else
    subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
        14587, 14586};
end
varNamesToPlot = {'auc', 'latency', 'timeUnder', 'openTime', 'initVelocity', ...
    'closeTime', 'maxClosingVelocity', 'maxOpeningVelocity', 'blinkRate'};

% varIdxToUse = {[5 7 2 1 4 3 8 6 9],...
%     [1 3 4 5 7 6 8 2 9]};

varIdxToUse = {[2 7 5 1 3 4 6 8 9],...
    [2 7 5 1 3 4 6 8 9]};


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
for pp = 1:length(varNamesToPlot)

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
        ii = find(strcmp(varNamesToPlot{pp},allVarNames));
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
        slopes(ss, pp) = fitObj.Coefficients.Estimate(2);

        % Get the y value at median x and it will be our offset
        offsets(ss,pp) = fitObj.Coefficients.Estimate(2)*median(x)+fitObj.Coefficients.Estimate(1);

        % subject parameter data session 1
        y = sessOne.(allVarNames{ii});
        goodPoints = ~isnan(y);
        x = log10(sessOne.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weightsSessOne);
        slopesSessOne(ss, pp) = fitObj.Coefficients.Estimate(2);

        % Get the y value at median x and it will be our offset
        offsetsSessOne(ss,pp) = fitObj.Coefficients.Estimate(2)*median(x)+fitObj.Coefficients.Estimate(1);

        % subject parameter data session 2
        y = sessTwo.(allVarNames{ii});
        goodPoints = ~isnan(y);
        x = log10(sessTwo.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weightsSessTwo);
        slopesSessTwo(ss, pp) = fitObj.Coefficients.Estimate(2);

        % Get the y value at median x and it will be our offset
        offsetsSessTwo(ss, pp) = fitObj.Coefficients.Estimate(2)*median(x)+fitObj.Coefficients.Estimate(1);
    
    end

end

%% Make test-retest analysis.

% First we do 2 PCA analyses on the slope and offset session aggregate
% variables. We then project session 1 and session 2 vectors onto the PC1
% and PC2 vectors and correlate their new scores.

% Save measurements and names in cells
allMeasures = {slopes, offsets};
allMeasuresSessOne = {slopesSessOne,offsetsSessOne};
allMeasuresSessTwo = {slopesSessTwo,offsetsSessTwo};
allMeasureNames = {'slopes', 'offsets'};

% Loop over the set [slopes, offsets, (slopes offsets)], and for each data set perform an NMF
% decomposition. Save the nmfResults struct and plot some diagnostics

figure

for ii = 1:2

    % Standardize the measures
    standardized = (allMeasures{ii}-mean(allMeasures{ii}))./std(allMeasures{ii});
    sessionOneStandard = (allMeasuresSessOne{ii}-mean(allMeasuresSessOne{ii}))./std(allMeasuresSessOne{ii});
    sessionTwoStandard = (allMeasuresSessTwo{ii}-mean(allMeasuresSessTwo{ii}))./std(allMeasuresSessTwo{ii});

    % Re-order the variables for this measure to improve plotting
    standardized = standardized(:,varIdxToUse{ii});
    sessionOneStandard = sessionOneStandard(:,varIdxToUse{ii});
    sessionTwoStandard = sessionTwoStandard(:,varIdxToUse{ii});
    varNames = strcat(varNamesToPlot(varIdxToUse{ii}),allMeasureNames{ii});

    % Initialize the coeff variable
    coeff = zeros(9,3);

    % Handle sign reversal
    if ii==1
        idxToFlip = contains(varNames,{'latency'});
        coeff(1:3,1) = 1/sqrt(3); % velocity
        coeff(4:6,2) = 1/sqrt(3); % depth
    end

    if ii==2
        idxToFlip = contains(varNames,{'latency'});
        coeff(2:6,1) = 1/sqrt(5); % velocity
        coeff(7:8,2) = 1/sqrt(2); % depth
    end

    standardized(:,idxToFlip) = -standardized(:,idxToFlip);
    sessionOneStandard(:,idxToFlip) = -sessionOneStandard(:,idxToFlip);
    sessionTwoStandard(:,idxToFlip) = -sessionTwoStandard(:,idxToFlip);
    varNames(idxToFlip) = strcat(varNames(idxToFlip),'_Neg');

    fSparse = length(coeff(:));
    for nn=1:1000
        [~,iterH] = nnmf(standardized,3,'H0',coeff');
        penaltyMat = iterH;
        penaltyMat(iterH>0.5) = -(iterH(iterH>0.5)-1);
        thisSparse = norm(penaltyMat(:));
        if thisSparse < fSparse && ~any(all(iterH'==0))
            fSparse = thisSparse;
            H = iterH;
        end
    end
    coeff = H';

    % If we are dealing with slopes, switch the order of the NMF dimensions
    % to make plotting a bit easier
    if ii==1
        coeff=coeff(:,[2 1 3]);
    end

    dColors = {'b','g'};

    subplot(3,2,ii)
    imagesc(corr(standardized,'Type','Kendall'))
    xticks(1:9);
    yticks(1:9);
    xticklabels(varNames);
    yticklabels(varNames);
    title([allMeasureNames{ii} ' correlation matrix'])
    for dd=1:2
        xx = find(coeff(:,dd)>0.35,1,"first");
        hw = find(coeff(:,dd)>0.35,1,"last")-xx+1;
        rectangle('Position',[xx-0.5 xx-0.5 hw hw],'EdgeColor',dColors{dd},'LineWidth',2);
    end

    sessionBothProjected = standardized*coeff(:,1:2);
    sessionOneProjected = sessionOneStandard*coeff(:,1:2);
    sessionTwoProjected = sessionTwoStandard*coeff(:,1:2);

    % Store the results
    results.(allMeasureNames{ii}) = sessionBothProjected;

    zRange = ceil(max([sessionOneProjected(:); sessionTwoProjected(:)]));
    for dd=1:2
        subplot(3,2,ii+dd*2)
        mdl = fitlm(sessionOneProjected(:,dd), sessionTwoProjected(:,dd),'RobustOpts','on');
        plot(mdl, 'Marker', 'o', 'MarkerEdgeColor','none', 'MarkerFaceColor',dColors{dd})
        legend off
        xlabel(sprintf('Session 1 projected onto d%d',dd))
        ylabel(sprintf('Session 2 projected onto d%d',dd))
        axis equal
        xlim([-zRange zRange])
        ylim([-zRange zRange])
        axis square
        title([allMeasureNames{ii} sprintf(' d%d',dd)]);
    end

end

    warning(warnState);
