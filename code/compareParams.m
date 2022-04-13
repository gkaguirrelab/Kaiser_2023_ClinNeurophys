%% compareParams
% This script loads a blink data set into a MATLAB table variable. When
% run, it will aggregate data for a given subject for the given parameters 
% across sessions. It will create a matrix of slopes for each parameter.
%%

% load file path
dataPath = fileparts(fileparts(mfilename('fullpath')));
spreadsheet ='UPenn Ipsi Summary_25ms_02062022.csv';

% choose subject and parameters
subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    14587, 14586};
varNamesToPlot = {'aucI', 'latencyI', 'timeUnderI', 'openTimeI', 'initVelocityI', ...
    'closeTimeI', 'maxClosingVelocityI', 'maxOpeningVelocityI', 'excursionI', 'closuresI'};

% create MATLAB table variable
T = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames = T.Properties.VariableNames;
slopes = zeros(length(subList), length(varNamesToPlot));
slopesSessOne = zeros(length(subList), length(varNamesToPlot));
slopesSessTwo = zeros(length(subList), length(varNamesToPlot));

%% create slopes matrix containing the slope values for each var and subject
for vv = 1:length(varNamesToPlot)

    for ss = 1:length(subList)

        % find scans for desired subject
        scans = T(ismember(T.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);
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
        
        % subject parameter data session 1
        y = sessOne.(allVarNames{ii});
        goodPoints = ~isnan(y);
        x = log10(scans.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weightsSessOne);
        slopesSessTwo(ss, vv) = fitObj.Coefficients.Estimate(2);
        
        % subject parameter data session 2
        y = sessTwo.(allVarNames{ii});
        goodPoints = ~isnan(y);
        x = log10(scans.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weightsSessTwo);
        slopesSessOne(ss, vv) = fitObj.Coefficients.Estimate(2);

    end

end

comp = corr(slopes);

%% Perform a PCA analysis after standardizing the slope measures
X = (slopes-mean(slopes))./std(slopes);
[coeff,score,latent,tsquared,explained,mu] = pca(X);
figure
plot(explained)
xlabel('component'); ylabel('percent variance explained');
figure
biplot(coeff(:,1:2),'scores',score(:,1:2),'varLabels',varNamesToPlot)
axis equal
figure
scatter3(score(:,1),score(:,2),score(:,3),'or')
axis equal
xlabel('component 1'); ylabel('component 2'); zlabel('component 3');

%% PC1 correlation between session 1 and session 2

weights = abs(coeff(:,1)');

% get session one scores
sO = (slopesSessOne-mean(slopesSessOne))./std(slopesSessOne);
[coeff,score,latent,tsquared,explained,mu] = pca(sO, 'VariableWeights', weights);
xx = score(:,1);

% get session two scores
sT = (slopesSessTwo-mean(slopesSessTwo))./std(slopesSessTwo);
[coeff,score,latent,tsquared,explained,mu] = pca(sT, 'VariableWeights', weights);
yy = score(:,1);

figure();
scatter(xx,yy)
fitObj = fitlm(xx,yy,'RobustOpts', 'on');
hold on
plot(xx,fitObj.Fitted,'-r')
rsquare = fitObj.Rsquared.Ordinary;
if rsquare > 1 || rsquare < 0
    rsquare = nan;
end
title(['PC1 scores by session - ' sprintf(' R^2=%2.2f',rsquare)])
xlabel(['Session one'])
ylabel(['Session two'])
axis square
ylim(xlim)

%% compare scores to slopes

figure();
plotNum = 1;

for vv = 1:length(varNamesToPlot)
    subplot(2, length(varNamesToPlot), plotNum)
    scatter(xx,slopesSessOne(:,vv))
    fitObj = fitlm(xx,slopesSessOne(:,vv),'RobustOpts', 'on');
    hold on
    plot(xx,fitObj.Fitted,'-r')
    rsquare = fitObj.Rsquared.Ordinary;
    if rsquare > 1 || rsquare < 0
        rsquare = nan;
    end
    title(['PC1 scores and ' varNamesToPlot(vv) ' - ' sprintf(' R^2=%2.2f',rsquare)])
    xlabel(['Session one PC1 scores'])
    ylabel(['Session one slopes'])

    subplot(2, length(varNamesToPlot), plotNum + length(varNamesToPlot))
    scatter(yy,slopesSessTwo(:,vv))
    fitObj = fitlm(yy,slopesSessTwo(:,vv),'RobustOpts', 'on');
    hold on
    plot(yy,fitObj.Fitted,'-r')
    rsquare = fitObj.Rsquared.Ordinary;
    if rsquare > 1 || rsquare < 0
        rsquare = nan;
    end
    title(['PC1 scores and ' varNamesToPlot(vv) ' - ' sprintf(' R^2=%2.2f',rsquare)])
    xlabel(['Session two PC1 scores'])
    ylabel(['Session two slopes'])
    plotNum = plotNum + 1;

end