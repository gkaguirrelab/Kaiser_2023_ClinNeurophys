%% checkValuesAllTrials
% This script loads a blink data set into a MATLAB table variable. When
% run, it calculates the mean value for a blink feature across trials in an
% acquisition. It then plots each blink value, the means, and the residuals
% across puff intensities.
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
subList = {15512};
varNamesToPlot = {'auc'};

% create MATLAB table variable
T = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames = T.Properties.VariableNames;

%% calculate and plot mean resuidual values
for vv = 1:length(varNamesToPlot)
    
    figure();
    plotNum = 0;
    
    for ss = 1:length(subList)
        
        allbx = [];
        allby = [];
        plotNum = plotNum + 1;
        acqMeans = NaN(1,25);
        ii = find(strcmp(varNamesToPlot{vv},allVarNames));

        % find scans for desired subject
        scans = T(ismember(T.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);

        % get individual and mean values for each scan
        for zz = 1:25
            temp = scans(ismember(scans.scanNumber, zz+1),:);
            if ~isempty(temp)
                trialNums = (temp.stimIndex)'*zz;
                allbx = horzcat(allbx, trialNums);
                allby = horzcat(allby, (temp.(varNamesToPlot{vv}))');
                acqMeans(zz) = mean(temp.(varNamesToPlot{vv}), 'omitnan');
             end
        end
        
        % plot all puffs
        subplot(3,length(subList),plotNum);
        scatter(allbx,allby);
        title({['Subject ' num2str(subList{ss})]}, 'FontSize', 14)
        if plotNum == 1
            ylabel(['Individual puff ' varNamesToPlot{vv}], 'FontSize', 14)
            xlabel('Trial number across session', 'FontSize', 14)
        end
        
        % plot fits
        subplot(3,length(subList),plotNum+length(subList));
        fitObj = fitlm(allbx,allby,'RobustOpts', 'on');
        modelY = fitObj.Fitted;
        scatter(allbx,modelY, '_');
        if plotNum == 1
            ylabel([varNamesToPlot{vv} ' model fit value'], 'FontSize', 14)
        end
        
        % plot residuals
        residuals = allby - modelY';
        subplot(3,length(subList),plotNum+length(subList)*2);
        scatter(allbx,residuals);
        if plotNum == 1
            ylabel([varNamesToPlot{vv} ' residuals'], 'FontSize', 14)
        end
        
    end
        
end