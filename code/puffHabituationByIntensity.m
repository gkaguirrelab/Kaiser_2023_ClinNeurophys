%% puffHabituationByIntensity
% This script loads a blink data set into a MATLAB table variable. When
% run, it calculates the mean value for a blink feature across trials in an
% acquisition for a specific puff intensity. It then fits the mean 
% acquisition value across puff intensity.
% For a set of single trials within each, acquisition the modeled acquisition
% value is subtracted from the value of the measure for each trial to get
% the residuals. We then plot the mean residual value from each acquisition
% as a function of acquisition number.
% 
%       Scan PSI index (out of 26 scans, discarding scan 1):
%          3.5 PSI: [3 8 13 24 25]
%          7.5 PSI: [9 11 12 20 22]
%          15 PSI: [4 7 16 17 21]
%          30 PSI: [2 10 15 18 26]
%          60 PSI: [5 6 14 19 23]
%
%% set up parameters

% load file path
dataPath = fileparts(fileparts(mfilename('fullpath')));
spreadsheet ='Upenn_Ipsilateral Afiles_clean_full.csv';

% choose subject and parameters
subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    14590, 14589, 14588, 14587, 14586};
varNamesToPlot = {'auc', 'latency', 'timeUnder20', 'openTime', 'initialVelocity', ...
     'closeTime', 'maxClosingVelocity', 'maxOpeningVelocity', 'blinkRate'};
psi = 30;
if psi == 3.5
    scanNums = [3 8 13 24 25];
elseif psi == 7.5
    scanNums = [9 11 12 20 22];
elseif psi == 15
    scanNums = [4 7 16 17 21];
elseif psi == 30
    scanNums = [2 10 15 18 26];
else
    scaNums = [5 6 14 19 23];
end

% create MATLAB table variable
T = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames = T.Properties.VariableNames;

%% calculate and plot mean resuidual values
for vv = 1:length(varNamesToPlot)
    
    figure();
    plotNum = 0;
    
    for ss = 1:length(subList)
        
        plotNum = plotNum + 1;
        acqMeans = NaN(1,5);
        resMeansByAcq = NaN(1,5);
        ii = find(strcmp(varNamesToPlot{vv},allVarNames));

        % find scans for desired subject
        scans = T(ismember(T.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);

        % get mean values for each scan
        for zz = 1:5
            temp = scans(ismember(scans.scanNumber, scanNums(zz)),:);
            if ~isempty(temp)
               acqMeans(zz) = mean(temp.(varNamesToPlot{vv}), 'omitnan');
            end
        end
        
        values = zeros(1,5) + psi;
        fitObj = fitlm(values,acqMeans,'RobustOpts', 'on');
        modelY = fitObj.Fitted;

        % calculate residuals as a function of acquisition number
        for zz = 1:5            
            temp = scans(ismember(scans.scanNumber, scanNums(zz)),:);            
            if ~isempty(temp)
                residuals = temp.(varNamesToPlot{vv})' - modelY(zz);
                resMeansByAcq(zz) = mean(residuals,'omitnan');
            end            
        end
        
        % calculate residuals as a function of trial number
        resByTrial = NaN(25, 8);
        resMeansByTrial = NaN(1,8);
        for zz = 1:8
            temp = scans(ismember(scans.stimIndex, zz),:);
            if ~isempty(temp)
                for yy = 1:5
                    tt = temp(ismember(temp.scanNumber, scanNums(yy)),:);
                    if isempty(tt)
                        residual = NaN;
                    elseif length(tt.(varNamesToPlot{vv})) == 1
                        residual = tt.(varNamesToPlot{vv})(1) - modelY(yy);
                    else
                        res1 = tt.(varNamesToPlot{vv})(1) - modelY(yy);
                        res2 = tt.(varNamesToPlot{vv})(1) - modelY(yy);
                        residual = mean([res1 res2]);
                    end
                    resByTrial(yy,zz) = residual;
                end
            end
            resMeansByTrial(1,zz) = mean(resByTrial(:,zz), 'omitnan');
        end

        % plot mean residual as a function of acquisition number
        subplot(2,length(subList),plotNum);
        scatter(scanNums,resMeansByAcq);
        fitObj = fitlm(scanNums,resMeansByAcq,'RobustOpts', 'on');
        hold on
        plot(scanNums,fitObj.Fitted,'-r');
        rsquare = fitObj.Rsquared.Ordinary;
        if rsquare > 1 || rsquare < 0
            rsquare = nan;
        end
        title({['Subject ' num2str(subList{ss})] ' when psi = ' psi,sprintf(' R^2=%2.2f',rsquare)}, 'FontSize', 14)
        if plotNum ~= 1
            yticklabels("");
            xticklabels("");
            xticks([]);
            yticks([]);
        else
            ylabel([varNamesToPlot{vv} ' residual'], 'FontSize', 14)
            xlabel('Acquisition number', 'FontSize', 14)
        end
        
        % plot mean residual as a function of trial number
        subplot(2,length(subList),plotNum + length(subList));
        scatter((1:8),resMeansByTrial);
        fitObj = fitlm((1:8),resMeansByTrial,'RobustOpts', 'on');
        hold on
        plot((1:8),fitObj.Fitted,'-r');
        rsquare = fitObj.Rsquared.Ordinary;
        if rsquare > 1 || rsquare < 0
            rsquare = nan;
        end
        title({['Subject ' num2str(subList{ss})] ' when psi = ' psi,sprintf(' R^2=%2.2f',rsquare)}, 'FontSize', 14)
        if plotNum ~= 1
            yticklabels("");
            xticklabels("");
            xticks([]);
            yticks([]);
        else
            ylabel([varNamesToPlot{vv} ' residual'], 'FontSize', 14)
            xlabel('Trial number', 'FontSize', 14)
        end

    end
        
end
