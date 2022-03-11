%% checkValues
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
subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    14590, 14589, 14588, 14587, 14586};
varNamesToPlot = {'auc', 'latency', 'timeUnder20', 'openTime', 'initialVelocity', ...
     'closeTime', 'maxClosingVelocity', 'maxOpeningVelocity', 'blinkRate'};

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
        A = [];
        B = [];
        C = [];
        D = [];
        E = [];
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
                if ismember(zz+1, [3 8 13 24 25])
                    psi = 3.5;
                    A = horzcat(A, (temp.(varNamesToPlot{vv}))');
                    psiArray = zeros(1, length(temp.(varNamesToPlot{vv}))) + psi;
                    allbx = horzcat(allbx, psiArray);
                    allby = horzcat(allby, (temp.(varNamesToPlot{vv}))');
                elseif ismember(zz+1, [9 11 12 20 22])
                    psi = 7.5;
                    B = horzcat(B, (temp.(varNamesToPlot{vv}))');
                    psiArray = zeros(1, length(temp.(varNamesToPlot{vv}))) + psi;
                    allbx = horzcat(allbx, psiArray);
                    allby = horzcat(allby, (temp.(varNamesToPlot{vv}))');
                elseif ismember(zz+1, [4 7 16 17 21])
                    psi = 15;
                    C = horzcat(C, (temp.(varNamesToPlot{vv}))');
                    psiArray = zeros(1, length(temp.(varNamesToPlot{vv}))) + psi;
                    allbx = horzcat(allbx, psiArray);
                    allby = horzcat(allby, (temp.(varNamesToPlot{vv}))');
                elseif ismember(zz+1, [2 10 15 18 26])
                    psi = 30;
                    D = horzcat(D, (temp.(varNamesToPlot{vv}))');
                    psiArray = zeros(1, length(temp.(varNamesToPlot{vv}))) + psi;
                    allbx = horzcat(allbx, psiArray);
                    allby = horzcat(allby, (temp.(varNamesToPlot{vv}))');
                else
                    psi = 60;
                    E = horzcat(E, (temp.(varNamesToPlot{vv}))');
                    psiArray = zeros(1, length(temp.(varNamesToPlot{vv}))) + psi;
                    allbx = horzcat(allbx, psiArray);
                    allby = horzcat(allby, (temp.(varNamesToPlot{vv}))');
                end
                acqMeans(zz) = mean(temp.(varNamesToPlot{vv}), 'omitnan');
             end
        end
        
        scanpsi = [30 3.75 15 60 60 15 3.75 7.5 30 7.5 7.5 3.75 60 30 15 15 30 60 7.5 15 7.5 60 3.75 3.75 30];
        fitObj = fitlm(log10(scanpsi),acqMeans,'RobustOpts', 'on');
        modelY = fitObj.Fitted;
        
        % model fit x and ys by puff intensity
        modelsX = log10([3.5 7.5 15 30 60]);
        modelsY = [modelY(2) modelY(8) modelY(3) modelY(1) modelY(4)];
        
        % caclulate residuals
        resA = A - modelsY(1);
        resB = B - modelsY(2);
        resC = C - modelsY(3);
        resD = D - modelsY(4);
        resE = E - modelsY(5);
        
        % set up psi values
        valA = zeros(1,length(resA)) + 3.5;
        valB = zeros(1,length(resB)) + 7.5;
        valC = zeros(1,length(resC)) + 15;
        valD = zeros(1,length(resD)) + 30;
        valE = zeros(1,length(resE)) + 60;
        vals = horzcat(valA, valB, valC, valD, valE);
        
        % residuals x and ys
        resX = log10(vals);
        resY = horzcat(resA, resB, resC, resD, resE);

        % plot all puffs
        subplot(3,length(subList),plotNum);
        scatter(log10(allbx),allby);
        title({['Subject ' num2str(subList{ss})]}, 'FontSize', 14)
        if plotNum == 1
            ylabel(['Individual puff ' varNamesToPlot{vv}], 'FontSize', 14)
            xlabel('Puff intensity [log psi]', 'FontSize', 14)
        end
        
        % plot fits
        subplot(3,length(subList),plotNum+length(subList));
        scatter(modelsX,modelsY, '_');
        if plotNum == 1
            ylabel([varNamesToPlot{vv} ' model fit value'], 'FontSize', 14)
        end
        
        % plot residuals
        subplot(3,length(subList),plotNum+length(subList)*2);
        scatter(resX,resY);
        if plotNum == 1
            ylabel([varNamesToPlot{vv} ' residuals'], 'FontSize', 14)
        end
        
    end
        
end