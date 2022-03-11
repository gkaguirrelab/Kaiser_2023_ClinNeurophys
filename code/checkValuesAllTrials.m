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
session = 1;

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
        psi = [30 3.75 15 60 60 15 3.75 7.5 30 7.5 7.5 3.75 60 30 15 15 30 60 7.5 15 7.5 60 3.75 3.75 30];
        acqMeans = NaN(1,25);
        ii = find(strcmp(varNamesToPlot{vv},allVarNames));

        % find scans for desired subject
        scans = T(ismember(T.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);
        dates = unique(scans.scanDate);
        if session == 1
            scans = scans(ismember(scans.scanDate,dates(1,1)),:);
        else
            scans = scans(ismember(scans.scanDate,dates(2,1)),:);
        end

        % get individual and mean values for each scan
        for zz = 1:25
            temp = scans(ismember(scans.scanNumber, zz+1),:);
            if ~isempty(temp)
                trialNums = (temp.stimIndex)' + 8*(zz-1);
                allbx = horzcat(allbx, trialNums);
                allby = horzcat(allby, (temp.(varNamesToPlot{vv}))');
                acqMeans(zz) = mean(temp.(varNamesToPlot{vv}), 'omitnan');
             end
        end
        
        % get model values as a function of puff intensity
        fitObj = fitlm(log10(psi),acqMeans,'RobustOpts', 'on');
        modelFit = fitObj.Fitted;
        
        % calculate model fit points
        modelY = zeros(1,length(allbx));
        aTrials = horzcat((9:16),(49:56),(89:96),(177:192));
        bTrials = horzcat((57:64),(73:88),(145:152),(161:168));
        cTrials = horzcat((17:24),(41:48),(113:128),(153:160));
        dTrials = horzcat((1:8),(65:72),(105:112),(129:136),(193:200));
        eTrials = horzcat((25:40),(97:104),(137:144),(169:176));
        for i = 1:length(allbx)
            if ismember(allbx(i),aTrials)
                modelY(i) = modelFit(2);
            elseif ismember(allbx(i),bTrials)
                modelY(i) = modelFit(8);
            elseif ismember(allbx(i),cTrials)
                modelY(i) = modelFit(3);
            elseif ismember(allbx(i),dTrials)
                modelY(i) = modelFit(1);
            elseif ismember(allbx(i),eTrials)
                modelY(i) = modelFit(4);
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
        scatter(allbx,modelY, '_');
        if plotNum == 1
            ylabel([varNamesToPlot{vv} ' model fit value'], 'FontSize', 14)
        end
        
        % plot residuals
        residuals = allby - modelY;
        subplot(3,length(subList),plotNum+length(subList)*2);
        scatter(allbx,residuals);
        if plotNum == 1
            ylabel([varNamesToPlot{vv} ' residuals'], 'FontSize', 14)
        end
        
    end
        
end