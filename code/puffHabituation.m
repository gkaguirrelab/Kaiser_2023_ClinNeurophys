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
spreadsheet ='Upenn_Ipsilateral Afiles_clean_03032022.csv';

% choose subject and parameters
% all subjects = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
%    14587, 14586}
subList = {15512, 14596, 14594, 14592, 14592, 14591, 14590, 14589, 14587, 14586};
varNamesToPlot = {'latency'};

% create MATLAB table variable
T = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames = T.Properties.VariableNames;

%% calculate and plot mean resuidual values
for vv = 1:length(varNamesToPlot)
    
    figure();
    plotNum = 0;
    
    for ss = 1:length(subList)
        
        plotNum = plotNum + 1;
        acqMeans = NaN(1,25);
        resMeansByAcq = NaN(1,25);
        psi = [30 3.75 15 60 60 15 3.75 7.5 30 7.5 7.5 3.75 60 30 15 15 30 60 7.5 15 7.5 60 3.75 3.75 30];
        ii = find(strcmp(varNamesToPlot{vv},allVarNames));

        % find scans for desired subject
        scans = T(ismember(T.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);
        xScans = unique(scans.scanNumber);
        xPuffs = unique(scans.stimIndex);

        % get mean values for each scan
        for zz = 1:25
            temp = scans(ismember(scans.scanNumber, zz+1),:);
            if ~isempty(temp) 
               acqMeans(zz) = mean(temp.(varNamesToPlot{vv}), 'omitnan');
             end
        end

        fitObj = fitlm(log10(psi),acqMeans,'RobustOpts', 'on');
        modelY = fitObj.Fitted;

        % calculate residuals as a function of acquisition number
        for zz = 1:25            
            temp = scans(ismember(scans.scanNumber, zz+1),:);            
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
                for yy = 1:25
                    tt = temp(ismember(temp.scanNumber, yy+1),:);
                    if isempty(tt)
                        residual = NaN;
                    elseif length(tt.(varNamesToPlot{vv})) == 1
                        residual = tt.(varNamesToPlot{vv})(1) - modelY(yy);
                    else
                        residual = mean(tt.(varNamesToPlot{vv})(1) + tt.(varNamesToPlot{vv})(2) - 2*modelY(yy));
                    end
                    resByTrial(yy,zz) = residual;
                end
            end
            resMeansByTrial(1,zz) = mean(resByTrial(:,zz), 'omitnan');
        end

        % plot mean residual as a function of acquisition number
        subplot(2,length(subList),plotNum);
        scatter((1:25),resMeansByAcq);
        fitObj = fitlm((1:25),resMeansByAcq,'RobustOpts', 'on');
        hold on
        plot((1:25),fitObj.Fitted,'-r');
        rsquare = fitObj.Rsquared.Ordinary;
        if rsquare > 1 || rsquare < 0
            rsquare = nan;
        end
        title(['Subject ' num2str(subList{ss}) ' ' sprintf(' R^2=%2.2f',rsquare)], 'FontSize', 14)
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
        scatter((1:8),resMeansByTrial);
        fitObj = fitlm((1:8),resMeansByTrial,'RobustOpts', 'on');
        hold on
        plot((1:8),fitObj.Fitted,'-r');
        rsquare = fitObj.Rsquared.Ordinary;
        if rsquare > 1 || rsquare < 0
            rsquare = nan;
        end
        title(['Subject ' num2str(subList{ss}) ' ' sprintf(' R^2=%2.2f',rsquare)], 'FontSize', 14)
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