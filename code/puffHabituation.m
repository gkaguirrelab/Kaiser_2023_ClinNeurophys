%% puffHabituation
% This script loads a blink data set into a MATLAB table variable. When
% run, it calculates the mean value for a blink feature across trials in an
% acquisition. It then fits the mean acquisition value across puff intensity.
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
spreadsheet ='UPenn_AFiles_IpsiOnly_03022022.csv';

% choose subject and parameters
% all subjects = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
%    14587, 14586}
% only run subjects 149590, 14589, and 14588 for highest 3 PSI levels
subList = {14587};
varNamesToPlot = {'latency'};
highestOnly = false;

xFit = linspace(log10(3),log10(70),50);

% create MATLAB table variable
T = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames = T.Properties.VariableNames;

%% calculate and plot mean resuidual values
for vv = 1:length(varNamesToPlot)
    
    figure();
    plotNum = 0;
    
    for ss = 1:length(subList)
        
        plotNum = plotNum + 1;
        acqMeans = [];
        resMeansByAcq = [];
        psi = [];
        ii = find(strcmp(varNamesToPlot{vv},allVarNames));

        % find scans for desired subject
        scans = T(ismember(T.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);
        xScans = unique(scans.scanNumber);
        xPuffs = unique(scans.stimIndex);

        % get mean values for each scan
        for zz = 1:length(xScans)
            temp = scans(ismember(scans.scanNumber, zz+1),:);
            if ~isempty(temp) 
               acqMeans(end+1) = mean(temp.(varNamesToPlot{vv}), 'omitnan');
               if ismember(zz+1, [3 8 13 24 25])
                   psi(end+1) = 3.5;
               elseif ismember(zz+1, [9 11 12 20 22])
                   psi(end+1) = 7.5;
               elseif ismember(zz+1, [4 7 16 17 21])
                   psi(end+1) = 15;
               elseif ismember(zz+1, [2 10 15 18 26])
                   psi(end+1) = 30;
               else
                   psi(end+1) = 60;
               end
             end
        end

        fitObj = fitlm(log10(psi),acqMeans,'RobustOpts', 'on');
        modelY = fitObj.Fitted;

        % calculate residuals as a function of acquisition number
        for zz = 2:26            
            temp = scans(ismember(scans.scanNumber, zz),:);            
            if ~isempty(temp)
                residuals = temp.(varNamesToPlot{vv})' - modelY(zz-1);
                resMeansByAcq(end+1) = mean(residuals,'omitnan');
            end            
        end
        
        % calculate residuals as a function of trial number
        resByTrial = NaN(length(xScans), length(xPuffs));
        resMeansByTrial = NaN(1,length(xPuffs));
        for zz = 1:length(xPuffs)            
            temp = scans(ismember(scans.stimIndex, zz),:);
            if ~isempty(temp)
                for yy = 1:length(xScans)
                    scan = xScans(yy);
                    tt = temp(ismember(temp.scanNumber, scan),:);
                    if isempty(tt)
                        residual = NaN;
                    else
                        residual = tt.(varNamesToPlot{vv})(1) - modelY(yy);
                    end
                    resByTrial(yy,zz) = residual;
                end
            end
            resMeansByTrial(1,zz) = mean(resByTrial(:,zz), 'omitnan');
        end

        % plot mean residual as a function of acquisition number
        subplot(2,length(subList),plotNum);
        scatter(xScans,resMeansByAcq);
        fitObj = fitlm(xScans,resMeansByAcq,'RobustOpts', 'on');
        hold on
        plot(xScans,fitObj.Fitted,'-r');
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
            ylabel(['Session one ' varNamesToPlot{vv} ' residual'], 'FontSize', 14)
            xlabel('Acquisition number', 'FontSize', 14)
        end
        
        % plot mean residual as a function of trial number
        subplot(2,length(subList),plotNum + length(subList));
        scatter(xPuffs,resMeansByTrial);
        fitObj = fitlm(xPuffs,resMeansByTrial,'RobustOpts', 'on');
        hold on
        plot(xPuffs,fitObj.Fitted,'-r');
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
            ylabel(['Session one ' varNamesToPlot{vv} ' residual'], 'FontSize', 14)
            xlabel('Trial number', 'FontSize', 14)
        end

    end
        
end