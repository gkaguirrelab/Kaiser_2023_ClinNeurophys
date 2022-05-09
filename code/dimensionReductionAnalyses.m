%% dimensionReductionAnalyses
% This script loads a blink data set into a MATLAB table variable. When
% run, it will aggregate data for a given subject for the given parameters
% across sessions. The routine then conducts an analysis of the
% dimensionality of the measures. A NMF (favored for sparsity) is
% performed, and the test / re-test reliability of the session data is then
% examined after projection to these dimensions.

% Housekeeping
clear; close all

% Some constants for the analysis
nDimensions = 3;
groupCellThresh = 0.25;
offsetX = log10(15); % The PSI value at which the offset is calculated

% We use a constant state of the rng given that the NMF is stochastic
rng('default')

% Store the warning state and turn off some that will arise
warnState = warning();
warning('off','stats:nnmf:LowRank');
warning('off','stats:statrobustfit:IterationLimit');

% Load file path
dataPath = fileparts(fileparts(mfilename('fullpath')));
spreadsheet ='UPENN Summary with IPSI Responses_02072022_SquintCheck.csv';

% Choose subject and parameters. "Highest only" is the entire set of
% subjects
highestOnly = false;
if highestOnly
    %         subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    %             14590, 14589, 14588, 14587, 14586};
    subList = {15512, 15507, 15506, 15505, 14595, 14594, 14593, 14592, 14591, ...
        14590, 14589, 14588, 14587, 14586};
else
    subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
        14587, 14586};
end

subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    14590, 14589, 14588, 14587, 14586};


nSubs = length(subList);


% The names of the blink features
varNamesToPlot = {'auc', 'latency', 'timeUnder', 'openTime', 'initVelocity', ...
    'closeTime', 'maxClosingVelocity', 'maxOpeningVelocity', 'blinkRate'};

% These vector are used to re-order the blink features in plotting to more
% easily show their grouping into dimensions.
varIdxToUse = {[8 2 7 5 4 1 3 6],...
    [8 2 7 5 4 1 3 6]};

% The set of intended PSI values
intendedPSI = [3.5,7.5,15,30,60];

% create MATLAB table variable
T = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames = T.Properties.VariableNames;
slopes = zeros(length(subList), length(varNamesToPlot));
slopesSessOne = zeros(length(subList), length(varNamesToPlot));
slopesSessTwo = zeros(length(subList), length(varNamesToPlot));
offsets = zeros(length(subList), length(varNamesToPlot));
offsetsSessOne = zeros(length(subList), length(varNamesToPlot));
offsetsSessTwo = zeros(length(subList), length(varNamesToPlot));
meanResiduals = nan(length(subList),length(varNamesToPlot),length(intendedPSI));

%% create slopes matrix containing the slope values for each var and subject
for pp = 1:length(varNamesToPlot)

    for ss = 1:length(subList)

        % find scans for desired subject
        scans = T(ismember(T.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);
        scans = scans(ismember(scans.notSquint,'TRUE'),:);
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
        offsets(ss,pp) = fitObj.Coefficients.Estimate(2)*offsetX+fitObj.Coefficients.Estimate(1);

        % Store the mean residual value by intended PSI
        for kk=1:length(intendedPSI)
            thisIdx = find(scans.intendedPSI==intendedPSI(kk));
            if ~isempty(thisIdx)
                meanResiduals(ss,pp,kk) = mean(fitObj.Residuals.Raw(thisIdx));
            end
        end

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
        offsetsSessOne(ss,pp) = fitObj.Coefficients.Estimate(2)*offsetX+fitObj.Coefficients.Estimate(1);

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
        offsetsSessTwo(ss, pp) = fitObj.Coefficients.Estimate(2)*offsetX+fitObj.Coefficients.Estimate(1);

    end

end

% Save measurements and names in cells
allMeasures = {slopes, offsets};
allMeasuresSessOne = {slopesSessOne,offsetsSessOne};
allMeasuresSessTwo = {slopesSessTwo,offsetsSessTwo};
allMeasureNames = {'slopes', 'offsets'};


%% Residuals
% Make a plot of the residuals for each measure across pressure
figure
for pp=1:length(varNamesToPlot)
    subplot(3,3,pp);
    plot(repmat(log10(intendedPSI),nSubs,1),squeeze(meanResiduals(:,pp,:)),'.k');
    hold on
    plot(log10(intendedPSI),nanmean(squeeze(meanResiduals(:,pp,:))),'or','MarkerFaceColor','r');
    plot(log10(intendedPSI),nanmean(squeeze(meanResiduals(:,pp,:))),'-r');
    title(varNamesToPlot{pp})
end




%% Dimensionality and reproducibility
% Loop over the set [slopes, offsets, (slopes offsets)], and for each data set perform an NMF
% decomposition. Save the nmfResults struct and plot some diagnostics

% Prepare a figure
figure

% Loop over the two parameters of the model fit (slope and offset)
for ii = 1:2

    % Z-score standardize the measures
    standardized = (allMeasures{ii}-mean(allMeasures{ii}))./std(allMeasures{ii});
    sessionOneStandard = (allMeasuresSessOne{ii}-mean(allMeasuresSessOne{ii}))./std(allMeasuresSessOne{ii});
    sessionTwoStandard = (allMeasuresSessTwo{ii}-mean(allMeasuresSessTwo{ii}))./std(allMeasuresSessTwo{ii});

    % Re-order the variables to improve plotting
    standardized = standardized(:,varIdxToUse{ii});
    sessionOneStandard = sessionOneStandard(:,varIdxToUse{ii});
    sessionTwoStandard = sessionTwoStandard(:,varIdxToUse{ii});
    varNames = strcat(varNamesToPlot(varIdxToUse{ii}),allMeasureNames{ii});

    % Initialize the coeff variable
    coeff = zeros(length(varNames),nDimensions);

    % Depending upon the parameter (slope or offset), we do two things
    % here:
    % 1) identify which blink features should be sign inverted. This is
    % done because the predicted relationship between stimulus intensity /
    % sensitivity and response is inverted for these. E.g., the more
    % sensitive someone is, the shorter their latency value should be.
    % Close time also has this property.
    % 2) define a coefficient matrix to initialize the NMF solution

    if ii==1
        idxToFlip = contains(varNames,{'latency','closeTime'});
        coeff(1:4,1) = 1/sqrt(4); % velocity
        coeff(5:6,2) = 1/sqrt(2); % depth
        coeff(7:8,2) = 1/sqrt(2); % depth
    end

    if ii==2
        idxToFlip = contains(varNames,{'latency','closeTime'});
        coeff(2:7,1) = 1/sqrt(5); % overall blink response (speed and size)
        coeff(1,2) = 1/sqrt(1); % overall blink response (speed and size)
    end

    % Apply the idxToFlip (sign inversion)
    standardized(:,idxToFlip) = -standardized(:,idxToFlip);
    sessionOneStandard(:,idxToFlip) = -sessionOneStandard(:,idxToFlip);
    sessionTwoStandard(:,idxToFlip) = -sessionTwoStandard(:,idxToFlip);
    varNames(idxToFlip) = strcat(varNames(idxToFlip),'_Neg');

    % Obtain an NMF solution. The NMF solution is stochastic. We search
    % across many such solutions, and favor the solution that groups the
    % cells into the categories defined in coeff above
    fVal = 1e6;
    for nn=1:1000
        [~,iterH] = nnmf(standardized,nDimensions,'H0',coeff');
        iterFval = norm(coeff-round(iterH)');
        if iterFval < fVal && ~any(all(iterH'==0))
            fVal = iterFval;
            H = iterH;
        end
    end
    coeff = H';

    % Show the correlation matrices of the blink features
    dColors = {'b','g'};
    subplot(4,2,ii)
    corrMap = corr(standardized,'Type','Kendall');
    corrMap(1:(length(varNames)+1):end)=nan;
    imAlpha=ones(size(corrMap));
    imAlpha(isnan(corrMap))=0;
    imagesc(corrMap,'AlphaData',imAlpha);
    set(gca,'color',0*[1 1 1]);
    title([allMeasureNames{ii} ' correlation matrix'])
    colormap(redblue);
    cMax = ceil(max(abs(corrMap(:)))*10)/10;
    caxis([-cMax cMax]);
    axis square
    colorbar

    % Add some rectangles to indicate the cell groups into dimensions
    for dd=1:2
        xx = find(coeff(:,dd)>groupCellThresh,1,"first");
        hw = find(coeff(:,dd)>groupCellThresh,1,"last")-xx+1;
        rectangle('Position',[xx-0.5 xx-0.5 hw hw],'EdgeColor',dColors{dd},'LineWidth',2);
    end

    % Make a plot of the dimension weights
    subplot(4,2,2+ii)
    b=bar(coeff(:,1:2));
    ylabel('coefficient');
    b(1).FaceColor=dColors{1}; b(2).FaceColor=dColors{2};
    axis square

    % Project the separate sessions onto the first two dimensions of the
    % coefficients
    sessionBothProjected = standardized*coeff(:,1:2);
    sessionOneProjected = sessionOneStandard*coeff(:,1:2);
    sessionTwoProjected = sessionTwoStandard*coeff(:,1:2);

    % Store the results
    results.(allMeasureNames{ii}) = sessionBothProjected;

    % Find the range of projected values to use for setting plot limits
    zRange = ceil(max([sessionOneProjected(:); sessionTwoProjected(:)]));

    % Plot test / retest plots for the two parameters (slope, offset)
    for dd=1:2
        subplot(4,2,2+ii+dd*2)
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

% Restore the warning state
warning(warnState);
