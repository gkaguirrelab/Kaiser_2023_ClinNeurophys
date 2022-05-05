%% runPCAAnalysis
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
varNamesToPlot = {'aucI', 'latencyI', 'timeUnderI', 'openTimeI', 'initVelocityI', ...
     'closeTimeI', 'maxClosingVelocityI', 'maxOpeningVelocityI', 'blinkRate'};

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

% Do the PCA analysis, save a pcaResults struct and plot some diagnostics 
for ii = 1:2
    % Do PCA with the aggregate measurements
    standardized = (allMeasures{1,ii}-mean(allMeasures{1,ii}))./std(allMeasures{1,ii});
    [coeff,score,latent,tsquared,explained,mu] = pca(standardized);
    
    % Save PCA results to a structure
    pcaResults.(['variable' allMeasureNames{1,ii}]).('coeff') = coeff;
    pcaResults.(['variable' allMeasureNames{1,ii}]).('score') = score;
    
    % Plot the PCA related stuff 
    set(0,'CurrentFigure',figure1)
    subplot(2, 3, plotCounter)
    plot(explained)
    xlabel('component'); ylabel('percent variance explained');
    title(allMeasureNames{1,ii})
    subplot(2, 3, plotCounter+1)
    biplot(coeff(:,1:2),'scores',score(:,1:2),'varLabels',varNamesToPlot)
    title(allMeasureNames{1,ii})
    axis equal
    subplot(2, 3, plotCounter+2)
    scatter3(score(:,1),score(:,2),score(:,3),'or')
    axis square
    title(allMeasureNames{1,ii})
    xlabel('component 1'); ylabel('component 2'); zlabel('component 3');
    
    % Project the individual sessions to PCA 
    sessionOneStandard = (allMeasures{ii+1,1}-mean(allMeasures{ii+1,1}))./std(allMeasures{ii+1,1});
    sessionTwoStandard = (allMeasures{ii+1,2}-mean(allMeasures{ii+1,2}))./std(allMeasures{ii+1,2});
    sessionOneProjected = sessionOneStandard*coeff(:,1:2);
    sessionTwoProjected = sessionTwoStandard*coeff(:,1:2);    
    combined = [sessionOneProjected sessionTwoProjected];
    combinedSessions = [combinedSessions; combined];
    
    set(0,'CurrentFigure',figure2) 
    subplot(2,2,plotCounter2)
    mdl = fitlm(sessionOneProjected(:,1), sessionTwoProjected(:,1));
    plot(mdl, 'Marker', 'o', 'MarkerEdgeColor','b', 'MarkerFaceColor','b')
    legend off
    xlabel('Session 1 projected onto PC1')
    ylabel('Session 2 projected onto PC1')
    title(allMeasureNames{ii})
    axis equal
    xlim([-4 4])
    ylim([-4 4])
    axis square

    subplot(2,2,plotCounter2+1)
    mdl = fitlm(sessionOneProjected(:,2), sessionTwoProjected(:,2));
    plot(mdl, 'Marker', 'o', 'MarkerEdgeColor','b', 'MarkerFaceColor','b')
    legend off
    xlabel('Session 1 projected onto PC2')
    ylabel('Session 2 projected onto PC2')
    title(allMeasureNames{ii})    
    axis equal
    xlim([-4 4])
    ylim([-4 4])
    axis square    
    
    plotCounter = plotCounter + 3;
    plotCounter2 = plotCounter2 + 2;    
end

% Combine session one and two projected to make one test/retest plot
figure('Renderer', 'painters', 'Position', [164 71 1401 891])
subplot(1, 2, 1)
mdl = fitlm(combinedSessions(:,1), combinedSessions(:,3));
plot(mdl, 'Marker', 'o', 'MarkerEdgeColor','b', 'MarkerFaceColor','b');
legend off
title('Test/retest reliability for slopes and offset combined (PC1)')
xlabel('Session 1')
ylabel('Session 2')
axis equal
xlim([-4 4])
ylim([-4 4])
axis square 
subplot(1, 2, 2)
mdl = fitlm(combinedSessions(:,2), combinedSessions(:,4));
plot(mdl, 'Marker', 'o', 'MarkerEdgeColor','b', 'MarkerFaceColor','b')
legend off
title('Test/retest reliability for slopes and offset combined (PC2)')
xlabel('Session 1')
ylabel('Session 2')
axis equal
xlim([-4 4])
ylim([-4 4])
axis square 

%% PC1 correlation of slopes and offset scores
   
pScores = pcaResults.variableSlopes.score(:,1);
oScores = pcaResults.variableOffsets.score(:,1);

co = corrcoef(pScores, oScores);
corrs = co(1,2)

bootstat = sort(bootstrp(1000,@corr,pScores,oScores));
CIL = bootstat(25)
CIH = bootstat(975)

%% Do a PCA correlation on the combined slopes and offsets
slopesAndOffsets = [slopes offsets];
slopesAndOffsets = (slopesAndOffsets-mean(slopesAndOffsets))./std(slopesAndOffsets);
[slopeOffsetCoeff,slopeOffsetScore,latent,tsquared,explained,mu] = pca(slopesAndOffsets);

newVarNamesToPlot = {'Slope_aucI', 'Slope_latencyI', 'Slope_timeUnderI', 'Slope_openTimeI', 'Slope_initVelocityI', ...
     'Slope_closeTimeI', 'Slope_maxClosingVelocityI', 'Slope_maxOpeningVelocityI', 'Slope_blinkRate', ...
     'Offset_aucI', 'Offset_latencyI', 'Offset_timeUnderI', 'Offset_openTimeI', 'Offset_initVelocityI', ...
     'Offset_closeTimeI', 'Offset_maxClosingVelocityI', 'Offset_maxOpeningVelocityI', 'Offset_blinkRate'};

% Plot the results
figure('Renderer', 'painters', 'Position', [164 71 1401 891])
subplot(1, 3, 1)
plot(explained)
xlabel('component'); ylabel('percent variance explained');
title('Variance explained for slope and offset combined')
subplot(1, 3, 2)
biplot(slopeOffsetCoeff(:,1:2),'scores',slopeOffsetScore(:,1:2),'varLabels',newVarNamesToPlot)
title('Loadings for slope and offset combined')
axis equal
subplot(1, 3, 3)
scatter3(slopeOffsetScore(:,1),slopeOffsetScore(:,2),slopeOffsetScore(:,3),'or')
axis square
title('Scores for variance explained')
xlabel('component 1'); ylabel('component 2'); zlabel('component 3');

% %% PC1 correlation between session 1 and session 2
% 
% weights = abs(coeff(:,1)');
% 
% % get session one scores
% sO = (slopesSessOne-mean(slopesSessOne))./std(slopesSessOne);
% [coeff,score,latent,tsquared,explained,mu] = pca(sO, 'VariableWeights', weights);
% xx = score(:,1);
% 
% % get session two scores
% sT = (slopesSessTwo-mean(slopesSessTwo))./std(slopesSessTwo);
% [coeff,score,latent,tsquared,explained,mu] = pca(sT, 'VariableWeights', weights);
% yy = score(:,1);
% 
% figure();
% scatter(xx,yy)
% fitObj = fitlm(xx,yy,'RobustOpts', 'on');
% hold on
% plot(xx,fitObj.Fitted,'-r')
% rsquare = fitObj.Rsquared.Ordinary;
% if rsquare > 1 || rsquare < 0
%     rsquare = nan;
% end
% title(['PC1 scores by session - ' sprintf(' R^2=%2.2f',rsquare)])
% xlabel(['Session one'])
% ylabel(['Session two'])
% axis square
% ylim(xlim)
% 
% %% compare scores to slopes
% 
% figure();
% plotNum = 1;
% 
% for vv = 1:length(varNamesToPlot)
%     subplot(2, length(varNamesToPlot), plotNum)
%     scatter(xx,slopesSessOne(:,vv))
%     fitObj = fitlm(xx,slopesSessOne(:,vv),'RobustOpts', 'on');
%     hold on
%     plot(xx,fitObj.Fitted,'-r')
%     rsquare = fitObj.Rsquared.Ordinary;
%     if rsquare > 1 || rsquare < 0
%         rsquare = nan;
%     end
%     title(['PC1 scores and ' varNamesToPlot(vv) ' - ' sprintf(' R^2=%2.2f',rsquare)])
%     xlabel(['Session one PC1 scores'])
%     ylabel(['Session one slopes'])
% 
%     subplot(2, length(varNamesToPlot), plotNum + length(varNamesToPlot))
%     scatter(yy,slopesSessTwo(:,vv))
%     fitObj = fitlm(yy,slopesSessTwo(:,vv),'RobustOpts', 'on');
%     hold on
%     plot(yy,fitObj.Fitted,'-r')
%     rsquare = fitObj.Rsquared.Ordinary;
%     if rsquare > 1 || rsquare < 0
%         rsquare = nan;
%     end
%     title(['PC1 scores and ' varNamesToPlot(vv) ' - ' sprintf(' R^2=%2.2f',rsquare)])
%     xlabel(['Session two PC1 scores'])
%     ylabel(['Session two slopes'])
%     plotNum = plotNum + 1;
% 
% end